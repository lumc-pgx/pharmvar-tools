import argparse
from itertools import combinations
from multiprocessing import Pool
import sys

from algebra import Relation
from algebra.relations.supremal_based import compare, find_supremal, spanning_variant
from algebra.utils import fasta_sequence
from algebra.variants import parse_hgvs, patch

from api import get_alleles, get_variants, get_version
from config import get_gene


# Global state for workers
ctx_worker_alleles = {}
ctx_worker_reference = ""


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def worker(args):
    lhs, rhs = args
    relation = compare(ctx_worker_reference, ctx_worker_alleles[lhs], ctx_worker_alleles[rhs])
    return lhs, rhs, relation


def main():
    parser = argparse.ArgumentParser(description="Calculate all relations of a gene")
    parser.add_argument("--gene", help="Gene to operate on", required=True)
    parser.add_argument("--version", help="Specify PharmVar version", default=get_version())
    parser.add_argument("--cores", type=int, help="Specify number of cores to run on", default=None)
    parser.add_argument("--disable-cache", help="Disable read and write from cache", action="store_true")
    args = parser.parse_args()

    try:
        gene_info = get_gene(args.gene)
    except KeyError:
        print(f"ERROR: Gene {args.gene} not in configuration!", file=sys.stderr)
        sys.exit(-1)

    ref_seq_id = gene_info["ng_ref_seq_id"]
    with open(f"data/{ref_seq_id}.fasta", encoding="utf-8") as file:
        reference = fasta_sequence(file.readlines())

    variants = get_variants(args.gene, ref_seq_id, args.version, not args.disable_cache)
    alleles = get_alleles(args.gene, ref_seq_id, args.version, not args.disable_cache)

    global ctx_worker_alleles
    global ctx_worker_reference
    ctx_worker_reference = reference

    for allele in alleles:
        try:
            allele_variants = [parse_hgvs(var["hgvs"], reference)[0] for var in allele["variants"]]
            observed = patch(reference, allele_variants)
            spanning = spanning_variant(reference, observed, allele_variants)
            supremal = find_supremal(reference, spanning)
            ctx_worker_alleles[allele["name"]] = supremal
        except ValueError as e:
            eprint(f"ERROR: allele {allele['name']} - {e}")

    for variant in variants:
        try:
            allele = parse_hgvs(variant["hgvs"], reference)
            supremal = find_supremal(reference, allele[0])
            ctx_worker_alleles[f"variant_{variant['id']}"] = supremal
        except ValueError as e:
            eprint(f"ERROR: variant {variant['hgvs']} - {e}")

    with Pool(args.cores) as pool:
        relations = pool.map(worker, combinations(ctx_worker_alleles, 2))
        for lhs, rhs, relation in relations:
            if relation != Relation.DISJOINT:
                print(lhs, rhs, relation.value)


if __name__ == "__main__":
    main()
