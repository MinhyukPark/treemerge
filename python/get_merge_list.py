import click
import dendropy
import networkx
import numpy

import njmerge2

@click.option("--start-tree", type=click.Path(exists=True), help="Tree on the full set of taxa")
@click.option("--tree-files", type=click.Path(exists=True), nargs=-1, help="Constraint trees")
@click.option("--output-prefix", type=click.Path(), help="Prefix for output files --- pair-merged trees and mst")
def get_merge_list(start_tree, tree_files, output_prefix):
    print(list(tree_files))
    # mst = njmerge2.tree_to_mst(strefile, treefiles)
    # graph = networkx.Graph(mst)
    # trees = []
    # for e in graph.edges():
    #     i, j = e
    #     tifile = treefiles[i]
    #     tjfile = treefiles[j]


if __name__ == "__main__":
    get_merge_list()
