from multiprocessing import Pool
import ntpath
import subprocess

import click
import dendropy
import numpy



def run_starting_trees_wrapper(args):
    return run_starting_trees(*args)


def run_starting_trees(merger_choice, guide_choice, starting_tree, current_starting_tree, pair_fasta):
    QUICKSCRIPTS = "./QuickScripts/"
    if(guide_choice == "induced"):
        with open(current_starting_tree + ".induce.out", "w") as stdout_f:
            with open(current_starting_tree + ".induce.err", "w") as stderr_f:
                subprocess.call(["/usr/bin/time", "-v", "python", QUICKSCRIPTS + "induce_tree.py", "--input-tree", starting_tree, "--sequence-file", pair_fasta, "--output-file", current_starting_tree], stdout=stdout_f, stderr=stderr_f)
    elif(guide_choice == "computed"):
        with open(current_starting_tree + ".iqtree.out", "w") as stdout_f:
            with open(current_starting_tree + ".iqtree.err", "w") as stderr_f:
                subprocess.call(["/usr/bin/time", "-v", "iqtree2", "-s", pair_fasta, "-nt", "AUTO", "-ntmax", "2", "-seed", "10101", "-m", "GTR+G", "-pre", current_starting_tree], stdout=stdout_f, stderr=stderr_f)
        with open(current_starting_tree + ".collapse.out", "w") as stdout_f:
            with open(current_starting_tree + ".collapse.err", "w") as stderr_f:
                subprocess.call(["/usr/bin/time", "-v", "python", QUICKSCRIPTS + "collapse_root.py", "--input-tree", current_starting_tree + ".treefile", "--output-file", current_starting_tree, "--hide-prefix"], stdout=stdout_f, stderr=stderr_f)

    if(merger_choice in ["constrained-inc", "njmerge"]):
        node_distance_matrix = current_starting_tree[:-len(".tree")] + "_node.mat"
        node_distance_matrix_taxlist = node_distance_matrix + "_taxlist"
        with open(current_starting_tree + ".tree_to_dist.out", "w") as stdout_f:
            with open(current_starting_tree + ".tree_to_dist.err", "w") as stderr_f:
                subprocess.call(["/usr/bin/time", "-v", "python", QUICKSCRIPTS + "tree_to_dist.py", "--input-filename", current_starting_tree, "--output-path", node_distance_matrix, "--dist-type", "node"], stdout=stdout_f, stderr=stderr_f)
        with open(current_starting_tree + ".mat_to_taxlist.out", "w") as stdout_f:
            with open(current_starting_tree + ".dist_to_taxlist.err", "w") as stderr_f:
                subprocess.call(["/usr/bin/time", "-v", "python", QUICKSCRIPTS + "mldist_to_treemerge.py", "--input-filename", node_distance_matrix, "--output-path", node_distance_matrix_taxlist], stdout=stdout_f, stderr=stderr_f)


@click.command()
@click.option("--starting-tree", type=click.Path(exists=True), help="Tree on the full set of taxa")
@click.option("--files-needed", type=click.Path(exists=True), help="Required pair-merged treefiles list")
@click.option("--merger-choice", type=click.Choice(["njmerge", "gtm", "constrained-inc"], case_sensitive=False))
@click.option("--guide-choice", type=click.Choice(["induced", "computed"], case_sensitive=False))
@click.option("--output-prefix", type=click.Path(), help="Prefix for output pair-merged trees")
def setup_merger(starting_tree, files_needed, merger_choice, guide_choice, output_prefix):
    QUICKSCRIPTS = "./QuickScripts/"
    starting_trees_args_arr = []
    with open(files_needed, "r") as f:
        for current_line in f:
            pair_tree_filename,pair_fasta,first_tree,second_tree = [current_filename.strip() for current_filename in current_line.split(" ")]
            pair_tree_dendropy = None
            current_starting_tree = pair_tree_filename[:-len(".tre")] + "-starting.tree"
            starting_trees_args_arr.append((merger_choice, guide_choice, starting_tree, current_starting_tree, pair_fasta))

    pool = Pool(processes=2)
    pool.map(run_starting_trees_wrapper, starting_trees_args_arr)

if __name__ == "__main__":
    setup_merger()
