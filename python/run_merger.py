from multiprocessing import Pool
import ntpath
import subprocess
import sys

import click
import dendropy
import numpy

QUICKSCRIPTS = "/opt/QuickScripts/"

try:
    import njmergepair
except ImportError:
    pass


def run_raxml_wrapper(args):
    return run_raxml(*args)

def run_merger_wrapper(args):
    return run_merger(*args)

def run_merger(merger_choice, first_tree, second_tree, pair_tree_filename):
    current_starting_tree = pair_tree_filename[:-len(".tre")] + "-starting.tree"
    node_distance_matrix = current_starting_tree[:-len(".tree")] + "_node.mat"
    node_distance_matrix_taxlist = node_distance_matrix + "_taxlist"
    merger_suffix = None
    if(merger_choice == "gtm"):
        merger_suffix = ".gtm"
    elif(merger_choice == "constrained-inc"):
        merger_suffix = ".cinc"
    elif(merger_choice == "njmerge"):
        merger_suffix = ".njmerge"

    with open(pair_tree_filename + merger_suffix + ".out", "w") as stdout_f:
        with open(pair_tree_filename + merger_suffix + ".err", "w") as stderr_f:
            if(merger_choice == "gtm"):
                subprocess.call(["/usr/bin/time", "-v", "python", "/opt/GTM/gtm.py", "-s", current_starting_tree, "-o", pair_tree_filename + merger_suffix + ".raw", "-t", first_tree, second_tree], stdout=stdout_f, stderr=stderr_f)
                with open(pair_tree_filename + merger_suffix + ".resolve.out", "w") as resolve_stdout_f:
                    with open(pair_tree_filename + merger_suffix + ".resolve.err", "w") as resolve_stderr_f:
                        subprocess.call(["/usr/bin/time", "-v", "python", QUICKSCRIPTS + "resolve_tree.py", "--input-tree", pair_tree_filename + merger_suffix + ".raw", "--output-file", pair_tree_filename + merger_suffix, "--hide-prefix"], stdout=resolve_stdout_f, stderr=resolve_stderr_f)
            elif(merger_choice in ["constrained-inc", "njmerge"]):
                node_distance_matrix = current_starting_tree[:-len(".tree")] + "_node.mat"
                node_distance_matrix_taxlist = node_distance_matrix + "_taxlist"

                if(merger_choice == "constrained-inc"):
                    first_tree_dup = pair_tree_filename + "." + first_tree.replace("/", "") + ".dup"
                    second_tree_dup = pair_tree_filename + "." + second_tree.replace("/", "") + ".dup"
                    current_starting_tree_dup = pair_tree_filename + "." + current_starting_tree.replace("/", "") + ".dup"
                    subprocess.call(["cp", first_tree, first_tree_dup])
                    subprocess.call(["cp", second_tree, second_tree_dup])
                    subprocess.call(["cp", current_starting_tree, current_starting_tree_dup])
                    subprocess.call(["/usr/bin/time", "-v", "constraint_inc", "-i", node_distance_matrix, "-o", pair_tree_filename + merger_suffix,  "-q", "subtree", "-g", current_starting_tree_dup, "-t", first_tree_dup, second_tree_dup], stdout=stdout_f, stderr=stderr_f)
                elif(merger_choice == "njmerge"):
                    subprocess.call(["/usr/bin/time", "-v", "python", "/u/sciteam/minhyuk2/git_repos/treemerge/python/njmergepair_wrapper.py", "--starting-tree", current_starting_tree, "--first-tree", first_tree, "--second-tree", second_tree, "--node-distance-matrix", node_distance_matrix, "--output-prefix", pair_tree_filename + merger_suffix], stdout=stdout_f, stderr=stderr_f)
            # subprocess.call(["python", QUICKSCRIPTS + "resolve_tree.py", "--input-tree", pair_tree_filename + merger_suffix, "--output-file", pair_tree_filename + merger_suffix + ".resolved", "--hide-prefix"])
            pair_tree_dendropy = dendropy.Tree.get(path=pair_tree_filename + merger_suffix, schema="newick")
            pair_tree_dendropy.is_rooted = False
            pair_tree_dendropy.collapse_basal_bifurcation(set_as_unrooted_tree=True)
            pair_tree_dendropy.update_bipartitions()
            pair_tree_dendropy.write(path=pair_tree_filename + merger_suffix + ".clean", schema="newick", suppress_rooting=True)


def run_raxml(merger_choice, pair_fasta, pair_tree_filename):
    merger_suffix = None
    if(merger_choice == "gtm"):
        merger_suffix = ".gtm"
    elif(merger_choice == "constrained-inc"):
        merger_suffix = ".cinc"
    elif(merger_choice == "njmerge"):
        merger_suffix = ".njmerge"

    # with open(pair_tree_filename + merger_suffix + ".trimal.out", "w") as stdout_f:
        # with open(pair_tree_filename + merger_suffix + ".trimal.err", "w") as stderr_f:
            # subprocess.call(["/usr/bin/time", "-v", "trimal", "-in", pair_fasta, "-out", pair_fasta + ".trimmed", "-gt", "0.05"], stdout=stdout_f, stderr=stderr_f)
    with open(pair_tree_filename + merger_suffix + ".raxml.out", "w") as stdout_f:
        with open(pair_tree_filename + merger_suffix + ".raxml.err", "w") as stderr_f:
            # subprocess.call(["/usr/bin/time", "-v", "raxmlng", "--evaluate", "--msa", pair_fasta + ".trimmed", "--threads", "2", "--model", "GTR+G", "--tree", pair_tree_filename + merger_suffix + ".clean", "--prefix", pair_tree_filename], stdout=stdout_f, stderr=stderr_f)
            subprocess.call(["/usr/bin/time", "-v", "raxmlng", "--evaluate", "--msa", pair_fasta, "--threads", "2", "--model", "GTR+G", "--tree", pair_tree_filename + merger_suffix + ".clean", "--prefix", pair_tree_filename], stdout=stdout_f, stderr=stderr_f)
            subprocess.call(["cp", pair_tree_filename + ".raxml.bestTree", pair_tree_filename])



@click.command()
@click.option("--starting-tree", type=click.Path(exists=True), help="Tree on the full set of taxa")
@click.option("--files-needed", type=click.Path(exists=True), help="Required pair-merged treefiles list")
@click.option("--merger-choice", type=click.Choice(["njmerge", "gtm", "constrained-inc"], case_sensitive=False))
@click.option("--guide-choice", type=click.Choice(["induced", "computed"], case_sensitive=False))
@click.option("--output-prefix", type=click.Path(), help="Prefix for output pair-merged trees")
def run_merger_and_raxml(starting_tree, files_needed, merger_choice, guide_choice, output_prefix):
    merger_args_arr = []
    raxml_args_arr = []
    with open(files_needed, "r") as f:
        for current_line in f:
            pair_tree_filename,pair_fasta,first_tree,second_tree = [current_filename.strip() for current_filename in current_line.split(" ")]
            pair_tree_dendropy = None
            merger_args_arr.append((merger_choice, first_tree, second_tree, pair_tree_filename))
            raxml_args_arr.append((merger_choice, pair_fasta, pair_tree_filename))

    pool = Pool(processes=30)
    pool.map(run_merger_wrapper, merger_args_arr)
    pool.close()

    pool = Pool(processes=15)
    pool.map(run_raxml_wrapper, raxml_args_arr)
    pool.close()

if __name__ == "__main__":
    run_merger_and_raxml()
