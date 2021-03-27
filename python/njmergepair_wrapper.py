import ntpath
import subprocess
import sys

import click
import dendropy
import numpy

try:
    import njmergepair
except ImportError:
    pass

@click.command()
@click.option("--starting-tree", type=click.Path(exists=True), help="Tree on the full set of taxa")
@click.option("--first-tree", type=click.Path(exists=True), help="First Tree")
@click.option("--second-tree", type=click.Path(exists=True), help="Second Tree")
@click.option("--node-distance-matrix", type=click.Path(exists=True), help="Node Distance Matrix")
@click.option("--output-prefix", type=click.Path(), help="Prefix for output pair-merged trees")
def njmergepair_wrapper(starting_tree, first_tree, second_tree, node_distance_matrix, output_prefix):
    QUICKSCRIPTS = "/projects/tallis/minhyuk2/git_repos/QuickScripts/"
    first_tree_dendropy = dendropy.Tree.get(path=first_tree, schema="newick")
    second_tree_dendropy = dendropy.Tree.get(path=second_tree, schema="newick")
    node_distance_matrix_taxlist = node_distance_matrix + "_taxlist"
    _,pair_tree_dendropy = njmergepair.run(node_distance_matrix, node_distance_matrix_taxlist, first_tree_dendropy, second_tree_dendropy)

    pair_tree_dendropy.resolve_polytomies()
    pair_tree_dendropy.update_bipartitions()
    pair_tree_dendropy.write(path=output_prefix, schema="newick", suppress_rooting=True)

if __name__ == "__main__":
    njmergepair_wrapper()
