import argparse
from itertools import combinations
import sys

from algebra import Relation
import networkx as nx

from api import get_version
import config


def read_relation(line):
    lhs, rhs, predicate, *_ = line.split()
    return lhs, rhs, predicate


def read_relations(file):
    relations = []
    for line in file.readlines():
        lhs, rhs, predicate = read_relation(line)
        relations.append((lhs, rhs, predicate))
    return relations


def write_edge(lhs, rhs, predicate):
    if predicate == Relation.IS_CONTAINED.value:
        return f'"{lhs}" -> "{rhs}";'
    if predicate == Relation.CONTAINS.value:
        return f'"{rhs}" -> "{lhs}";'

    rhs, lhs = sorted([lhs, rhs])
    if predicate == Relation.EQUIVALENT.value:
        return f'"{lhs}" -> "{rhs}" [arrowsize=0, color="black:invis:black"];'
    if predicate == Relation.OVERLAP.value:
        return f'"{lhs}" -> "{rhs}" [arrowsize=0, style=dashed];'
    if predicate != Relation.DISJOINT.value:
        raise ValueError(f"unknown relation: {predicate}")


def write_node(label, attributes):
    def quote(value):
        return f'"{value}"'

    return f'"{label}" [{",".join([key + "=" + quote(value) for key, value in attributes.items()])}];'


def write_dot(relations, nodes=None, file=sys.stdout):
    print("digraph {", file=file)
    if not nodes:
        nodes = {}
    print("\n".join([write_node(key, value) for key, value in nodes.items()]), file=file)
    print("\n".join([write_edge(*relation) for relation in relations]), file=file)
    print("}", file=file)


def build_graphs(relations):
    equivalent = nx.Graph()
    containment = nx.DiGraph()
    overlap = nx.Graph()
    for relation in relations:
        lhs, rhs, predicate = relation
        if predicate == Relation.EQUIVALENT.value:
            equivalent.add_edge(lhs, rhs)
        elif predicate == Relation.IS_CONTAINED.value:
            containment.add_edge(lhs, rhs)
        elif predicate == Relation.CONTAINS.value:
            containment.add_edge(rhs, lhs)
        elif predicate == Relation.OVERLAP.value:
            overlap.add_edge(lhs, rhs)
        elif predicate != Relation.DISJOINT.value:
            raise ValueError(f"unknown relation: {predicate}")

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
    for edge in overlap.edges():
        lhs, rhs = edge
        if (lhs not in ancestors or rhs not in ancestors or
                ancestors[lhs].isdisjoint(ancestors[rhs])):
            selected.add_edge(lhs, rhs)
    return selected


def most_specific_overlap(containment, overlap):
    to_remove = set()
    for node in overlap.nodes():
        for pair in combinations(overlap.edges(node), 2):
            (_, lhs), (_, rhs) = pair
            if lhs in containment.nodes() and rhs in containment.nodes():
                if lhs in nx.ancestors(containment, rhs):
                    to_remove.add((node, rhs))
                elif rhs in nx.ancestors(containment, lhs):
                    to_remove.add((node, lhs))

    overlap.remove_edges_from(to_remove)
    overlap.remove_nodes_from(set(nx.isolates(overlap)))
    return overlap


def export_relations(graph, predicate):
    relations = []
    for edge in graph.edges():
        relations.append((*edge, predicate))
    return relations


def select_context(equivalent, containment, overlap, context):
    nodes = set(context)
    for node in context:
        # TODO: no containment if there is already equivalence?
        if node in containment.nodes():
            nodes.update(list(nx.ancestors(containment, node)))
        # TODO: no overlap if there is already containment or equivalence?
        if node in overlap.nodes():
            nodes.update(list(overlap.neighbors(node)))

    context = set(nodes)
    for node in context:
        if node in equivalent.nodes():
            nodes.update(list(nx.node_connected_component(equivalent, node)))

    return equivalent.subgraph(nodes), containment.subgraph(nodes), overlap.subgraph(nodes)


def simplify(relations, context=None):
    equivalent, containment, overlap = build_graphs(relations)
    equivalent, containment, overlap = contract_equivalent(equivalent, containment, overlap)
    containment = nx.transitive_reduction(containment)
    overlap = overlap_without_common_ancestor(containment, overlap)
    overlap = most_specific_overlap(containment, overlap)

    if context:
        equivalent, containment, overlap = select_context(equivalent, containment, overlap, context)

    return equivalent, containment, overlap


def prepare4dot(equivalent, containment, overlap, nodes, context):
    for node in context:
        if node not in nodes:
            nodes[node] = {}

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
    parser.add_argument("--context", help="File with contextual nodes")
    parser.add_argument("--reference", help="Reference to operate on (default: %(default)s)", choices=["NG", "NC"], default="NG")
    parser.add_argument("--version", help="Specify PharmVar version", default=get_version())
    parser.add_argument("--disable-cache", help="Disable read and write from cache", action="store_true")

    args = parser.parse_args()

    try:
        gene_info = config.get_gene(args.gene)
    except KeyError:
        print(f"ERROR: Gene {args.gene} not in configuration!", file=sys.stderr)
        sys.exit(-1)

    if args.reference == "NG":
        ref_seq_id = gene_info["ng_ref_seq_id"]
    else:
        ref_seq_id = gene_info["nc_ref_seq_id"]

    nodes = config.get_nodes(args.gene, args.version, not args.disable_cache, ref_seq_id)

    context = set()
    if args.context:
        with open(args.context, encoding="utf-8") as file:
            context = set(file.read().splitlines())

    relations = read_relations(sys.stdin)
    write_dot(*prepare4dot(*simplify(relations, context), nodes, context))


if __name__ == "__main__":
    main()
