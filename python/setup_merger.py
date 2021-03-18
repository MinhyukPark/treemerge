import ntpath
import subprocess

import click
import dendropy
import numpy

@click.command()
@click.option("--starting-tree", type=click.Path(exists=True), help="Tree on the full set of taxa")
@click.option("--files-needed", type=click.Path(exists=True), help="Required pair-merged treefiles list")
@click.option("--merger-choice", type=click.Choice(["njmerge", "gtm", "constrained-inc"], case_sensitive=False))
@click.option("--guide-choice", type=click.Choice(["induced", "computed"], case_sensitive=False))
@click.option("--output-prefix", type=click.Path(), help="Prefix for output pair-merged trees")
def setup_merger(starting_tree, files_needed, merger_choice, guide_choice, output_prefix):
    QUICKSCRIPTS = "/projects/tallis/minhyuk2/git_repos/QuickScripts/"
    with open(files_needed, "r") as f:
        for current_line in f:
            pair_tree_filename,pair_fasta,first_tree,second_tree = [current_filename.strip() for current_filename in current_line.split(" ")]
            pair_tree_dendropy = None
            current_starting_tree = pair_tree_filename[:-len(".tre")] + "-starting.tree"

            if(guide_choice == "induced"):
                subprocess.call(["python", QUICKSCRIPTS + "induce_tree.py", "--input-tree", starting_tree, "--sequence-file", pair_fasta, "--output-file", current_starting_tree])
            elif(guide_choice == "computed"):
                subprocess.call(["iqtree2", "-s", pair_fasta, "-nt", "AUTO", "-ntmax", "2", "-seed", "10101", "-m", "GTR+G", "-pre", current_starting_tree])
                subprocess.call(["python", QUICKSCRIPTS + "resolve_tree.py", "--input-tree", current_starting_tree + ".treefile", "--output-file", current_starting_tree, "--hide-prefix"])
                # subprocess.call(["cp", current_starting_tree + ".treefile", current_starting_tree])

            if(merger_choice in ["constrained-inc", "njmerge"]):
                node_distance_matrix = current_starting_tree[:-len(".tree")] + "_node.mat"
                node_distance_matrix_taxlist = node_distance_matrix + "_taxlist"
                subprocess.call(["python", QUICKSCRIPTS + "tree_to_dist.py", "--input-filename", current_starting_tree, "--output-path", node_distance_matrix, "--dist-type", "node"])
                subprocess.call(["python", QUICKSCRIPTS + "mldist_to_treemerge.py", "--input-filename", node_distance_matrix, "--output-path", node_distance_matrix_taxlist])

if __name__ == "__main__":
    setup_merger()
