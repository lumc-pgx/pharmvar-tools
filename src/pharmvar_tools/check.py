import argparse
from itertools import combinations
from pathlib import Path
import requests
import sys

from algebra import Relation, are_equivalent, compare
from algebra.relations.sequence_based import are_equivalent as are_equivalent_sequence
from algebra.utils import fasta_sequence, vcf_variant
from algebra.variants import Variant, parse_hgvs, patch, reverse_complement
from ratelimiter import RateLimiter

from .api import get_alleles, get_variants, get_version
from .config import get_gene


NCBI_URI = "https://api.ncbi.nlm.nih.gov/variation/v0"
MUTALYZER_URI = "https://mutalyzer.nl/api"


def allele_from_variants(reference, variants):
    allele = []
    for variant in variants:
        allele.extend(parse_hgvs(variant["hgvs"], reference))
    return allele


def mutalyzer_hgvs(query):
    response = requests.get(f"{MUTALYZER_URI}/normalize/{query}")
    if response.status_code == 200:
        return response.json()["normalized_description"]
    return None


def check_alleles_vs_variants(alleles, variants, ref_seq_id):
    print(f"Checking variants with alleles for {ref_seq_id} ...")
    for allele in alleles:
        for variant in allele["variants"]:
            if variant not in variants:
                print(f"Variant {variant} from allele {allele['name']} not in variants.")

    for variant in variants:
        found = False
        for allele in alleles:
            if variant in allele["variants"]:
                found = True
        if not found:
            print(f"Variant {variant} not in alleles.")


def ncbi_hgvs(query):
    response = requests.get(f"{NCBI_URI}/hgvs/{query}/contextuals").json()

    position = int(response["data"]["spdis"][0]["position"])
    deleted = response["data"]["spdis"][0]["deleted_sequence"]
    inserted = response["data"]["spdis"][0]["inserted_sequence"]
    return Variant(position, position + len(deleted), inserted)


# NCBI Variation Services limited to 1 request/second.
@RateLimiter(max_calls=1, period=1)
def check_variants_with_ncbi(reference, ref_seq_id, variants):
    print(f"Checking variants with NCBI for {ref_seq_id} ...")
    for variant in variants:
        pharm_var = parse_hgvs(variant["hgvs"], reference)
        ncbi_var = [ncbi_hgvs(f"{ref_seq_id}:g.{variant['hgvs']}")]
        if pharm_var != ncbi_var and not are_equivalent(reference, pharm_var, ncbi_var):
            print(f"Non equivalent variants for {variant['hgvs']}: {pharm_var} vs {ncbi_var}")


def check_variants_with_mutalyzer(reference, ref_seq_id, variants):
    print(f"Checking variants with Mutalyzer for {ref_seq_id} ...")
    for variant in variants:
        pharm_var = parse_hgvs(variant["hgvs"], reference)
        mut_hgvs = mutalyzer_hgvs(f"{ref_seq_id}:g.{variant['hgvs']}")
        mut_var = parse_hgvs(mut_hgvs, reference)
        if pharm_var != mut_var and not are_equivalent(reference, pharm_var, mut_var):
            print(f"Non equivalent variants for {variant['hgvs']}: {pharm_var} vs {mut_var}")


def check_hgvs_allele_vs_variant_list(reference, ref_seq_id, alleles):
    print(f"Checking consistency between 'hgvs' and 'variants' entries in alleles for {ref_seq_id} ...")
    for allele in alleles:
        try:
            hgvs_var = parse_hgvs(allele["hgvs"], reference)
        except ValueError as error:
            print(f"Parsing of {allele['hgvs']} ({allele['name']}) failed ({str(error)})")
            continue

        variants = allele_from_variants(reference, allele["variants"])

        if hgvs_var != sorted(variants):
            print(f"HGVS variant {allele['hgvs']} mismatches with variant list ({allele['name']})")


def check_allele_variants(reference, ref_seq_id, alleles):
    print(f"Checking consistency of 'variants' entries of alleles for {ref_seq_id} ...")
    for allele in alleles:
        variants = set(allele_from_variants(reference, allele["variants"]))
        try:
            list(sorted(variants))
        except ValueError:
            print(f"{allele['name']} ({ref_seq_id}) unorderable")

            for lhs, rhs in combinations(variants, 2):
                relation = compare(reference, [lhs], [rhs])
                if relation != Relation.DISJOINT:
                    print(f"    {ref_seq_id}:g.{lhs.to_hgvs(reference)} and {ref_seq_id}:g.{rhs.to_hgvs(reference)} {relation.value} in {allele['name']}")


def check_allele_duplicates(reference, ref_seq_id, alleles):
    print(f"Checking for duplicates in 'variants' entries of alleles for {ref_seq_id} ...")
    for allele in alleles:
        variants = allele_from_variants(reference, allele["variants"])
        for duplicate in set([variant for variant in variants if variants.count(variant) > 1]):
            print(f"Duplicate {ref_seq_id}:g.{duplicate.to_hgvs(reference)} in {allele['name']}")


def multi_fasta(lines):
    name = None
    sequence = None
    for line in lines:
        if line.startswith(">"):
            if name is not None:
                yield name, sequence
            name, *_ = line[1:].split()
            sequence = ""
        else:
            sequence += line.strip()
    if name is not None:
        yield name, sequence


def check_hgvs_allele_vs_fasta(data_dir, reference, ref_seq_id, alleles, gene, version):
    print(f"Checking consistency between allele hgvs and fasta for {ref_seq_id} ...")
    fasta_alleles = {}
    with open(Path(data_dir, f"pharmvar-{version}", gene, f"{gene}.haplotypes.fasta"), encoding="utf-8") as file:
        for name, sequence in multi_fasta(file):
            fasta_alleles[name] = sequence
    for allele in alleles:
        try:
            hgvs_var = parse_hgvs(allele["hgvs"], reference)
        except ValueError as error:
            print(f"Parsing of {allele['hgvs']} ({allele['name']}) failed ({str(error)})")
            continue

        try:
            fasta_allele = fasta_alleles[allele["name"]]
        except KeyError:
            print(f"Fasta sequence not found for {allele['name']}")
            continue

        if not are_equivalent_sequence(reference, patch(reference, hgvs_var), fasta_allele):
            print(f"Non equivalent variants for {allele['hgvs']}: {hgvs_var} vs fasta ({allele['name']})")


def check_hgvs_allele_vs_vcf_ng(data_dir, gene, reference, ref_seq_id, alleles, version):
    print(f"Checking consistency between allele hgvs and NG vcf for {ref_seq_id} ...")
    for allele in alleles:
        try:
            hgvs_var = parse_hgvs(allele["hgvs"], reference)
        except ValueError as error:
            print(f"Parsing of {allele['hgvs']} ({allele['name']}) failed ({str(error)})")
            continue

        vcf_variants = []
        try:
            file = open(Path(data_dir, f"pharmvar-{version}", gene, "RefSeqGene", f"{allele['name'].replace('*', '_')}.vcf"), encoding="utf-8")
            for line in file:
                if not line.startswith("#"):
                    vcf_variants.append(vcf_variant(line))
            file.close()
        except FileNotFoundError:
            print(f"VCF file not found for {allele['name']}")
            continue

        if hgvs_var != vcf_variants and not are_equivalent(reference, hgvs_var, vcf_variants):
            print(f"Non equivalent variants for {allele['hgvs']}: {hgvs_var} vs {vcf_variants} ({allele['name']})")


def check_hgvs_allele_vs_vcf_nc(data_dir, gene, reference, ref_seq_id, alleles, version):
    print(f"Checking consistency between NC variants and NC vcf for {ref_seq_id} ...")
    for allele in alleles:
        try:
            hgvs_var = allele_from_variants(reference, allele["variants"])
        except ValueError as error:
            print(f"Parsing of {allele['name']} failed ({str(error)})")
            continue

        vcf_variants = []
        try:
            file = open(Path(data_dir, f"pharmvar-{version}", gene, "GRCh38", f"{allele['name'].replace('*', '_')}.vcf"), encoding="utf-8")
            for line in file:
                if not line.startswith("#"):
                    vcf_variants.append(vcf_variant(line))
            file.close()
        except FileNotFoundError:
            print(f"VCF file not found for {allele['name']}")
            continue

        try:
            if hgvs_var != vcf_variants and not are_equivalent(reference, hgvs_var, vcf_variants):
                print(f"Non equivalent variants for {allele['name']}: {hgvs_var} vs {vcf_variants}")
        except ValueError:
            # silently skip parsing/interpretation related problems checked elsewhere
            pass


def check_nc_vs_ng(nc_reference, ng_reference, nc_alleles, ng_alleles, mapping):
    print(f"Checking consistency between NC and NG ...")
    nc_mapped = nc_reference[mapping["start"] - 1:mapping["end"]]

    for nc_allele in nc_alleles:
        for ng_allele in ng_alleles:
            if nc_allele["name"] == ng_allele["name"]:
                try:
                    nc_variants = allele_from_variants(nc_reference, nc_allele["variants"])
                    ng_variants = allele_from_variants(ng_reference, ng_allele["variants"])

                    nc_observed = patch(nc_mapped, [Variant(variant.start - mapping["start"] + 1, variant.end - mapping["start"] + 1, variant.sequence) for variant in nc_variants])
                    if mapping["reverse"]:
                        nc_observed = reverse_complement(nc_observed)
                    ng_observed = patch(ng_reference, ng_variants)
                except ValueError:
                    # silently skip parsing/interpretation related problems checked elsewhere
                    break

                if nc_observed != ng_observed:
                    print(f"NC is not consistent with NG for {nc_allele['name']}")
                    print(nc_variants)
                    print(ng_variants)
                break


def main():
    parser = argparse.ArgumentParser(description="PharmVar data checker")
    parser.add_argument("--all", help="Perform all checks", action="store_true")
    parser.add_argument("--local", help="Perform all local checks", action="store_true")
    parser.add_argument("--alleles-vs-variants", help="Check allele vs. variant endpoints", action="store_true")
    parser.add_argument("--ncbi", help="Check variants against NCBI", action="store_true")
    parser.add_argument("--mutalyzer", help="Check variants against Mutalyzer", action="store_true")
    parser.add_argument("--variants", help="Check variants of alleles", action="store_true")
    parser.add_argument("--duplicates", help="Check for duplicate variants in alleles", action="store_true")
    parser.add_argument("--hgvs", help="Check allele hgvs entry vs. variant list", action="store_true")
    parser.add_argument("--fasta", help="Check hgvs entry of allele vs. fasta files", action="store_true")
    parser.add_argument("--vcf", help="Check hgvs entry of allele vs. vcf files", action="store_true")
    parser.add_argument("--nc-vs-ng", help="Check NC variants vs. NG variants", action="store_true")
    parser.add_argument("--gene", help="Gene to operate on", required=True)
    parser.add_argument("--disable-cache", help="Disable read and write from cache", action="store_true")
    parser.add_argument("--data-dir", help="Data directory", default="./data")
    parser.add_argument("--version", help="Specify PharmVar version")

    args = parser.parse_args()

    try:
        gene_info = get_gene(args.gene)
    except KeyError:
        print(f"ERROR: Gene {args.gene} not in configuration!")
        sys.exit(-1)

    nc_ref_seq_id = gene_info["nc_ref_seq_id"]
    ng_ref_seq_id = gene_info["ng_ref_seq_id"]


    print("Loading reference data ...")
    with open(Path(args.data_dir, f"{nc_ref_seq_id}.fasta"), encoding="utf-8") as file:
        nc_reference = fasta_sequence(file.readlines())

    with open(Path(args.data_dir, f"{ng_ref_seq_id}.fasta"), encoding="utf-8") as file:
        ng_reference = fasta_sequence(file.readlines())

    cache = not args.disable_cache
    print("Retrieving variant data ...")
    nc_variants = get_variants(args.data_dir, args.gene, nc_ref_seq_id, args.version, cache)
    ng_variants = get_variants(args.data_dir, args.gene, ng_ref_seq_id, args.version, cache)
    print("Retrieving allele data ...")
    nc_alleles = get_alleles(args.data_dir, args.gene, nc_ref_seq_id, args.version, cache)
    ng_alleles = get_alleles(args.data_dir, args.gene, ng_ref_seq_id, args.version, cache)

    if args.alleles_vs_variants or args.all or args.local:
        check_alleles_vs_variants(nc_alleles, nc_variants, nc_ref_seq_id)
        check_alleles_vs_variants(ng_alleles, ng_variants, ng_ref_seq_id)

    if args.ncbi or args.all:
        check_variants_with_ncbi(nc_reference, nc_ref_seq_id, nc_variants)
        check_variants_with_ncbi(ng_reference, ng_ref_seq_id, ng_variants)

    if args.mutalyzer or args.all:
        check_variants_with_mutalyzer(nc_reference, nc_ref_seq_id, nc_variants)
        check_variants_with_mutalyzer(ng_reference, ng_ref_seq_id, ng_variants)

    if args.variants or args.all or args.local:
        check_allele_variants(ng_reference, ng_ref_seq_id, ng_alleles)
        check_allele_variants(nc_reference, nc_ref_seq_id, nc_alleles)

    if args.duplicates or args.all or args.local:
        check_allele_duplicates(nc_reference, nc_ref_seq_id, nc_alleles)
        check_allele_duplicates(ng_reference, ng_ref_seq_id, ng_alleles)

    if not args.version:
        args.version = get_version()

    # Only for NG, as Allele["hgvs"] is always expressed as NG
    if args.hgvs or args.all or args.local:
        check_hgvs_allele_vs_variant_list(ng_reference, ng_ref_seq_id, ng_alleles)
    if args.fasta or args.all or args.local:
        check_hgvs_allele_vs_fasta(args.data_dir, ng_reference, ng_ref_seq_id, ng_alleles, args.gene, args.version)
    if args.vcf or args.all or args.local:
        check_hgvs_allele_vs_vcf_nc(args.data_dir, args.gene, nc_reference, nc_ref_seq_id, nc_alleles, args.version)
        check_hgvs_allele_vs_vcf_ng(args.data_dir, args.gene, ng_reference, ng_ref_seq_id, ng_alleles, args.version)

    if args.nc_vs_ng or args.all or args.local:
        check_nc_vs_ng(nc_reference, ng_reference, nc_alleles, ng_alleles, gene_info["nc_mapping"])


if __name__ == "__main__":
    main()
