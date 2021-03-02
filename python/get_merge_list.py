import ntpath

import click
import dendropy
import networkx
import numpy

import njmerge2

@click.command()
@click.option("--starting-tree", type=click.Path(exists=True), help="Tree on the full set of taxa")
@click.option("--output-prefix", type=click.Path(), help="Prefix for output files --- pair-merged trees and mst")
@click.argument("treefiles", type=click.Path(exists=True), nargs=-1)
def get_merge_list(starting_tree,  output_prefix, treefiles):
    """
    TREEFILES is the list of space separated filenames for the constraint trees
    """
    mst = njmerge2.tree_to_mst(starting_tree, treefiles)
    graph = networkx.Graph(mst)
    trees = []
    with open(output_prefix + "files_needed", "w") as f:
        for current_edge in graph.edges():
            i,j = current_edge
            tifile = treefiles[i]
            tjfile = treefiles[j]
            pair_tree_filename = ntpath.split(tifile)[1] + "+" + ntpath.split(tjfile)[1] + ".tree"
            f.write(pair_tree_filename)
            f.write("\n")
    networkx.write_gpickle(graph, output_prefix + "mst")


if __name__ == "__main__":
    get_merge_list()
