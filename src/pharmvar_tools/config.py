from importlib import resources
import yaml

from .api import get_alleles, get_variants


_COLORS = None
_GENES = None


def get_gene(gene):
    global _GENES
    if not _GENES:
        res = resources.files('pharmvar_tools') / 'config' / 'genes.yaml'
        with res.open(encoding="utf-8") as file:
            _GENES = yaml.safe_load(file)
    return _GENES[gene]


def function_to_color(function):
    global _COLORS
    if not _COLORS:
        res = resources.files('pharmvar_tools') / 'config' / 'colors.yaml'
        with res.open(encoding="utf-8") as file:
            _COLORS = yaml.safe_load(file)
    return _COLORS[function]


def impact_to_color(impact):
    if impact is None:
        return function_to_color("normal function")
    if any(["del" in impact, "dup" in impact, "ins" in impact]):
        return function_to_color("decreased function")
    if any(["fs" in impact, "splice" in impact, "X" in impact]):
        return function_to_color("no function")
    if impact:
        # FIXME: AA substitutions (BLOSOM?)
        return function_to_color("unknown function")
    # not assigned
    return function_to_color("function not assigned")


def get_nodes(data_dir, gene, version, cache, ref_seq_id):
    nodes = {}

    nodes.update({f"variant_{variant['id']}": {
        "color": impact_to_color(variant["impact"]),
        "label": variant["hgvs"],
        "shape": "oval",
    } for variant in get_variants(data_dir, gene, ref_seq_id, version, cache)})

    nodes.update({allele["name"]: {
        "color": function_to_color(allele["function"]),
        "label": allele["name"],
        "shape": "oval",
    } for allele in get_alleles(data_dir, gene, ref_seq_id, version, cache)})
    return nodes
