import argparse
from itertools import combinations
import sys

from algebra import Relation
import networkx as nx

from .api import get_version
from . import config


def read_relations(file=sys.stdin):
    return [line.split() for line in file.readlines()]


def dot_edge(lhs, rhs, value):
    if value == Relation.IS_CONTAINED.value:
        return f'"{lhs}" -> "{rhs}";'
    if value == Relation.CONTAINS.value:
        return f'"{rhs}" -> "{lhs}";'

    rhs, lhs = sorted([lhs, rhs])
    if value == Relation.EQUIVALENT.value:
        return f'"{lhs}" -> "{rhs}" [arrowsize=0, color="black:invis:black"];'
    if value == Relation.OVERLAP.value:
        return f'"{lhs}" -> "{rhs}" [arrowsize=0, style=dashed];'
    if value != Relation.DISJOINT.value:
        raise ValueError(f"unknown relation: {value}")


def quote(value):
    return f'"{value}"'


def dot_node(label, attributes):
    return f'"{label}" [{", ".join([key + "=" + quote(value) for key, value in attributes.items()])}];'


def write_dot(edges, nodes=None, sinks=None, file=sys.stdout):
    dot = create_dot(edges, nodes, sinks)
    print(dot, file=file)


def create_dot(edges, nodes=None, sinks=None):
    dot = "digraph {\n"
    if not nodes:
        nodes = {}
    dot += "\n".join([dot_node(key, value) for key, value in nodes.items()])
    dot += "\n"
    dot += "\n".join([dot_edge(*edge) for edge in edges])
    dot += "\n"
    if sinks:
        dot += f'{{rank="sink"; {"; ".join([quote(sink) for sink in sinks])}}}\n'
    dot += "}\n"

    return dot


def build_graphs(relations):
    equivalent = nx.Graph()
    containment = nx.DiGraph()
    overlap = nx.Graph()
    for relation in relations:
        lhs, rhs, value = relation
        if value == Relation.EQUIVALENT.value:
            equivalent.add_edge(lhs, rhs)
        elif value == Relation.IS_CONTAINED.value:
            containment.add_edge(lhs, rhs)
        elif value == Relation.CONTAINS.value:
            containment.add_edge(rhs, lhs)
        elif value == Relation.OVERLAP.value:
            overlap.add_edge(lhs, rhs)
        elif value != Relation.DISJOINT.value:
            raise ValueError(f"unknown relation: {value}")
    return equivalent, containment, overlap


def contract_equivalent(equivalent, containment, overlap):
    collapsed = nx.Graph()
    for component in nx.connected_components(equivalent):
        nodes = sorted(list(component))
        for node in nodes[1:]:
            if nodes[0] in containment.nodes() and node in containment.nodes():
                containment = nx.contracted_nodes(containment, nodes[0], node)
            if nodes[0] in overlap.nodes() and node in overlap.nodes():
                overlap = nx.contracted_nodes(overlap, nodes[0], node)
            collapsed.add_edge(nodes[0], node)
    return collapsed, containment, overlap


def overlap_without_common_ancestor(containment, overlap):
    ancestors = {}
    for node in containment.nodes():
        nodes = nx.ancestors(containment, node)
        if len(nodes) > 0:
            ancestors[node] = set(nodes)

    selected = nx.Graph()
    for lhs, rhs in overlap.edges():
        if (lhs not in ancestors or rhs not in ancestors or
                ancestors[lhs].isdisjoint(ancestors[rhs])):
            selected.add_edge(lhs, rhs)
    return selected


def most_specific_overlap(containment, overlap):
    to_remove = set()
    for node in overlap.nodes():
        for (_, lhs), (_, rhs) in combinations(overlap.edges(node), 2):
            if lhs in containment.nodes() and rhs in containment.nodes():
                if lhs in nx.ancestors(containment, rhs):
                    to_remove.add((node, rhs))
                elif rhs in nx.ancestors(containment, lhs):
                    to_remove.add((node, lhs))

    overlap.remove_edges_from(to_remove)
    overlap.remove_nodes_from(set(nx.isolates(overlap)))
    return overlap


def select_context(equivalent, containment, overlap, context):
    nodes = set(context)
    for node in context:
        if node in containment.nodes():
            nodes.update(list(nx.ancestors(containment, node)))
        if node in overlap.nodes():
            nodes.update(list(overlap.neighbors(node)))

    context = set(nodes)
    for node in context:
        if node in equivalent.nodes():
            nodes.update(list(nx.node_connected_component(equivalent, node)))

    return equivalent.subgraph(nodes), containment.subgraph(nodes), overlap.subgraph(nodes)


def simplify(equivalent, containment, overlap, context=None):
    equivalent, containment, overlap = contract_equivalent(equivalent, containment, overlap)
    containment = nx.transitive_reduction(containment)
    overlap = overlap_without_common_ancestor(containment, overlap)
    overlap = most_specific_overlap(containment, overlap)

    if context:
        equivalent, containment, overlap = select_context(equivalent, containment, overlap, context)

    return equivalent, containment, overlap


def export_relations(graph, value):
    relations = []
    for edge in graph.edges():
        relations.append((*edge, value))
    return relations


def prepare4export(equivalent, containment, overlap, nodes, context):
    for node in context:
        if node not in nodes:
            nodes[node] = {"shape": "box"}

    return (export_relations(equivalent, Relation.EQUIVALENT.value) +
            export_relations(containment, Relation.IS_CONTAINED.value) +
            export_relations(overlap, Relation.OVERLAP.value),
            {node: nodes[node] for node in nodes if node in
                list(equivalent.nodes()) +
                list(containment.nodes()) +
                list(overlap.nodes()) +
                list(context)})


def main():
    parser = argparse.ArgumentParser(description="Create Graphviz dot file of relations")
    parser.add_argument("--gene", help="Gene to operate on", required=True)
    parser.add_argument("--reference", help="Reference to operate on (default: %(default)s)", choices=["NG", "NC"], default="NG")
    parser.add_argument("--version", help="Specify PharmVar version")
    parser.add_argument("--disable-cache", help="Disable read and write from cache", action="store_true")
    parser.add_argument("--disable-simplify", help="Disable simplification of relations", action="store_true")
    parser.add_argument("--context", nargs='*', help="List of contextual nodes", default=[])
    parser.add_argument("--text", help="Plain text output", action="store_true")
    parser.add_argument("--data-dir", help="Data directory", default="./data")

    args = parser.parse_args()

    if not args.version:
        args.version = get_version()

    try:
        gene_info = config.get_gene(args.gene)
    except KeyError:
        print(f"ERROR: Gene {args.gene} not in configuration", file=sys.stderr)
        sys.exit(-1)

    if args.reference == "NG":
        ref_seq_id = gene_info["ng_ref_seq_id"]
    else:
        ref_seq_id = gene_info["nc_ref_seq_id"]

    config_nodes = config.get_nodes(args.data_dir, args.gene, args.version, not args.disable_cache, ref_seq_id)
    equivalent, containment, overlap = build_graphs(read_relations())
    if not args.disable_simplify:
        equivalent, containment, overlap = simplify(equivalent, containment, overlap, args.context)
    edges, dot_nodes = prepare4export(equivalent, containment, overlap, config_nodes, args.context)

    if args.text:
        for lhs, rhs, value in edges:
            print(f"{lhs} {rhs} {value}")
    else:
        write_dot(edges, dot_nodes, args.context)


if __name__ == "__main__":
    main()
