import argparse
from multiprocessing import Pool
import os
from pathlib import Path
import pickle
import sys

from algebra.lcs import LCSgraph
from algebra.utils import fasta_sequence
from algebra.variants import parse_hgvs

from .api import get_alleles, get_variants, get_version
from .config import get_gene


def read_lcsgraphs(path):
    with open(f"{path}", "rb") as file:
        return pickle.load(file)


def init_worker(*data):
    global worker_reference
    global worker_pv_alleles
    worker_reference, worker_pv_alleles = data


def worker(idx):
    try:
        variants = parse_hgvs(worker_pv_alleles[idx]["hgvs"] or "=", worker_reference)
    except ValueError:
        print(f"parsing failed for {worker_pv_alleles[idx]['name']}")
        return idx, None
    return idx, LCSgraph.from_variant(worker_reference, variants)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


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

    if args.lcsgraphs:
        graphs_file = args.lcsgraphs
    else:
        graphs_file = f"{args.data_dir}/pharmvar-{args.version}_{args.gene}_{ref_seq_id}_lcsgraphs.txt"
    if os.path.isfile(graphs_file):
        raise ValueError(f"LCS graphs file {graphs_file} exists already")

    pv_variants = get_variants(args.data_dir, args.gene, ref_seq_id, args.version, not args.disable_cache)
    pv_alleles = get_alleles(args.data_dir, args.gene, ref_seq_id, args.version, not args.disable_cache)

    alleles = {}
    print("Calculating LCS graphs for pharmvar alleles ...")
    with Pool(processes=args.cores, initializer=init_worker, initargs=(reference, pv_alleles)) as pool:
        result = pool.map(worker, range(len(pv_alleles)))
        for idx, graph in result:
            if graph is not None:
                alleles[pv_alleles[idx]["name"]] = graph

    print("Calculating LCS graphs for pharmvar variants ...")
    for allele in pv_variants:
        graph = LCSgraph.from_variant(reference, parse_hgvs(allele["hgvs"], reference))
        alleles[f"variant_{allele['id']}"] = graph

    with open(f"{graphs_file}", "wb") as file:
        pickle.dump(alleles, file)


if __name__ == "__main__":
    main()
