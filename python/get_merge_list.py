import click
import dendropy
import networkx
import numpy
import pickle

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
            f.write(tifile + "+" + tjfile)
            f.write("\n")
    pickle.write_gpickle(graph, output_prefix + "mst")


if __name__ == "__main__":
    get_merge_list()
