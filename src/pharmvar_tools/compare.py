import argparse
from itertools import combinations
from multiprocessing import Pool
import os
from pathlib import Path
import sys

from algebra import Relation
from algebra.relations.graph_based import compare
from algebra.utils import fasta_sequence
from algebra.variants import parse_spdi

from .api import get_alleles, get_variants, get_version
from .config import get_gene
from .lcsgraphs import read_lcsgraphs


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
    parser.add_argument("--lcsgraphs", help="File with LCS graphs to operate on")
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

    if args.lcsgraphs:
        lcsgraphs_file = args.lcsgraphs
    else:
        lcsgraphs_file = f"{args.data_dir}/pharmvar-{args.version}_{args.gene}_{ref_seq_id}_lcsgraphs.txt"
    if not os.path.isfile(lcsgraphs_file):
        raise ValueError(f"LCS graphs file {lcsgraphs_file} does not exist")

    lcsgraphs = read_lcsgraphs(lcsgraphs_file)

    alleles = {}
    for allele in pv_alleles:
        try:
            try:
                graph = lcsgraphs[allele["name"]]
            except KeyError:
                continue
            alleles[allele["name"]] = graph
        except ValueError as e:
            eprint(f"ERROR: allele {allele['name']} - {e}")

    for variant in pv_variants:
        try:
            try:
                graph = lcsgraphs[f"variant_{variant['id']}"]
            except KeyError:
                continue
            alleles[f"variant_{variant['id']}"] = graph
        except ValueError as e:
            eprint(f"ERROR: variant {variant['hgvs']} - {e}")

    with Pool(args.cores, initializer=init_worker, initargs=(reference, alleles)) as pool:
        relations = pool.map(worker, combinations(alleles, 2))
        for lhs, rhs, relation in relations:
            if relation != Relation.DISJOINT:
                print(lhs, rhs, relation.value)


if __name__ == "__main__":
    main()
