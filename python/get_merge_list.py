import ntpath
import subprocess

import click
import dendropy
import networkx
import numpy

import njmerge2
import treemerge

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
            pair_tree_filename = treemerge.name_treepair_file(output_prefix, tifile, tjfile)
            # pair_tree_filename = tifile + "+" + tjfile + ".tree"
            first_fasta = tifile[:-9] + ".out"
            second_fasta = tjfile[:-9] + ".out"
            merged_fasta_filename = pair_tree_filename + ".fasta"
            with open(merged_fasta_filename, "w") as fw:
                with open(first_fasta, "r") as ff:
                    for line in ff:
                        fw.write(line)
                with open(second_fasta, "r") as fs:
                    for line in fs:
                        fw.write(line)
            f.write(pair_tree_filename)
            f.write(" ")
            f.write(merged_fasta_filename)
            f.write(" ")
            f.write(tifile)
            f.write(" ")
            f.write(tjfile)
            f.write("\n")
    networkx.write_gpickle(graph, output_prefix + "mst")


if __name__ == "__main__":
    get_merge_list()
