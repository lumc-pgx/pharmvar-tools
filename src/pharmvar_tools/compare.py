import argparse
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path
import sys

from algebra import Relation
from algebra.relations.supremal_based import compare, find_supremal, spanning_variant
from algebra.utils import fasta_sequence
from algebra.variants import parse_hgvs, patch

from .api import get_alleles, get_variants, get_version
from .check import allele_from_variants
from .config import get_gene


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def init_worker(*data):
    global worker_reference
    global worker_alleles
    worker_reference, worker_alleles = data


def worker(args):
    lhs, rhs = args
    relation = compare(worker_reference, worker_alleles[lhs], worker_alleles[rhs])
    return lhs, rhs, relation


def main():
    parser = argparse.ArgumentParser(description="Calculate all relations of a gene")
    parser.add_argument("--gene", help="Gene to operate on", required=True)
    parser.add_argument("--reference", help="Reference to operate on (default: %(default)s)", choices=["NG", "NC"], default="NG")
    parser.add_argument("--version", help="Specify PharmVar version")
    parser.add_argument("--cores", type=int, help="Specify number of cores to run on", default=None)
    parser.add_argument("--data-dir", help="Data directory", default="./data")
    parser.add_argument("--disable-cache", help="Disable read and write from cache", action="store_true")
    args = parser.parse_args()

    if not args.version:
        args.version = get_version()

    try:
        gene_info = get_gene(args.gene)
    except KeyError:
        print(f"ERROR: Gene {args.gene} not in configuration!", file=sys.stderr)
        sys.exit(-1)

    if args.reference == "NG":
        ref_seq_id = gene_info["ng_ref_seq_id"]
    else:
        ref_seq_id = gene_info["nc_ref_seq_id"]

    with open(Path(args.data_dir, f"{ref_seq_id}.fasta"), encoding="utf-8") as file:
        reference = fasta_sequence(file.readlines())

    pv_variants = get_variants(args.data_dir, args.gene, ref_seq_id, args.version, not args.disable_cache)
    pv_alleles = get_alleles(args.data_dir, args.gene, ref_seq_id, args.version, not args.disable_cache)

    alleles = {}
    for allele in pv_alleles:
        try:
            allele_variants = allele_from_variants(reference, allele["variants"])
            observed = patch(reference, allele_variants)
            spanning = spanning_variant(reference, observed, allele_variants)
            supremal, *_ = find_supremal(reference, spanning)
            alleles[allele["name"]] = supremal
        except ValueError as e:
            eprint(f"ERROR: allele {allele['name']} - {e}")

    for variant in pv_variants:
        try:
            allele = parse_hgvs(variant["hgvs"], reference)
            supremal, *_ = find_supremal(reference, allele[0])
            alleles[f"variant_{variant['id']}"] = supremal
        except ValueError as e:
            eprint(f"ERROR: variant {variant['hgvs']} - {e}")

    with Pool(args.cores, initializer=init_worker, initargs=(reference, alleles)) as pool:
        relations = pool.map(worker, combinations(alleles, 2))
        for lhs, rhs, relation in relations:
            if relation != Relation.DISJOINT:
                print(lhs, rhs, relation.value)


if __name__ == "__main__":
    main()
