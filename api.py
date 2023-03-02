import json

import requests


PHARMVAR_URI = "https://www.pharmvar.org/api-service"


def get_version():
    url = "https://www.pharmvar.org/get-version"
    response = requests.get(url).json()
    return response["versionNumber"]


# FIXME in Python 3.9: .removeprefix()
def _remove_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text


def _to_variants(data, ref_seq_id=""):
    variants = []
    for variant in data:
        variants.append({
            "hgvs": _remove_prefix(variant["hgvs"], f"{ref_seq_id}:g."),
            "id": variant["variantId"],
            "impact": variant["impact"],
        })
    return variants


def _cache_requests(url, params, cache, path):
    if cache:
        try:
            with open(f"{path}", encoding="utf-8") as file:
                return json.load(file)
        except FileNotFoundError:
            pass

    response = requests.get(url, params=params).json()

    if cache:
        with open(f"{path}", "w", encoding="utf-8") as file:
            json.dump(response, file)

    return response


def get_alleles(gene, ref_seq_id, version, cache=False):
    response = _cache_requests(f"{PHARMVAR_URI}/genes/{gene}", {
                                   "reference-location-type": "Sequence Start",
                                   "reference-sequence": {ref_seq_id},
                               }, cache, f"data/pharmvar-{version}_{gene}_{ref_seq_id}_alleles.json")

    alleles = []
    for allele in response["alleles"]:
        alleles.append({
            "function": allele["function"],
            "hgvs": _remove_prefix(allele["hgvs"], f"{ref_seq_id}:g."),
            "name": allele["alleleName"],
            "variants": _to_variants(allele["variants"], ref_seq_id),
        })
    return alleles


def get_variants(gene, ref_seq_id, version, cache=False):
    response = _cache_requests(f"{PHARMVAR_URI}/variants/gene/{gene}", {
                                   "reference-location-type": "Sequence Start",
                                   "reference-sequence": {ref_seq_id},
                               }, cache, f"data/pharmvar-{version}_{gene}_{ref_seq_id}_variants.json")
    return _to_variants(response, ref_seq_id)
