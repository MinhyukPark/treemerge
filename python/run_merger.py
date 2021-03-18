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
@click.option("--files-needed", type=click.Path(exists=True), help="Required pair-merged treefiles list")
@click.option("--merger-choice", type=click.Choice(["njmerge", "gtm", "constrained-inc"], case_sensitive=False))
@click.option("--guide-choice", type=click.Choice(["induced", "computed"], case_sensitive=False))
@click.option("--output-prefix", type=click.Path(), help="Prefix for output pair-merged trees")
def run_merger_and_raxml(starting_tree, files_needed, merger_choice, guide_choice, output_prefix):
    QUICKSCRIPTS = "/projects/tallis/minhyuk2/git_repos/QuickScripts/"
    with open(files_needed, "r") as f:
        for current_line in f:
            pair_tree_filename,pair_fasta,first_tree,second_tree = [current_filename.strip() for current_filename in current_line.split(" ")]
            pair_tree_dendropy = None
            current_starting_tree = pair_tree_filename[:-len(".tre")] + "-starting.tree"
            merger_suffix = ""
            # first_tree_resolved = first_tree + ".resolved"
            # second_tree_resolved = second_tree + ".resolved"

            # subprocess.call(["python", QUICKSCRIPTS + "resolve_tree.py", "--input-tree", first_tree, "--output-file", first_tree_resolved, "--hide-prefix"])
            # subprocess.call(["python", QUICKSCRIPTS + "resolve_tree.py", "--input-tree", second_tree, "--output-file", second_tree_resolved, "--hide-prefix"])

            if(merger_choice == "gtm"):
                merger_suffix = ".gtm"
                subprocess.call(["python", "/opt/GTM/gtm.py", "-s", current_starting_tree, "-o", pair_tree_filename + merger_suffix, "-t", first_tree, second_tree])
                pair_tree_dendropy = dendropy.Tree.get(path=pair_tree_filename + merger_suffix, schema="newick")
            elif(merger_choice in ["constrained-inc", "njmerge"]):
                node_distance_matrix = current_starting_tree[:-len(".tree")] + "_node.mat"
                node_distance_matrix_taxlist = node_distance_matrix + "_taxlist"
                if(merger_choice == "constrained-inc"):
                    merger_suffix = ".cinc"
                    subprocess.call(["constraint_inc", "-i", node_distance_matrix, "-o", pair_tree_filename + merger_suffix,  "-q", "subtree", "-g", current_starting_tree, "-t", first_tree, second_tree])
                    # sys.stderr.write(" ".join(["constraint_inc", "-i", node_distance_matrix, "-o", pair_tree_filename + merger_suffix,  "-q", "subtree", "-g", current_starting_tree, "-t", first_tree, second_tree]))
                    pair_tree_dendropy = dendropy.Tree.get(path=pair_tree_filename + merger_suffix, schema="newick")
                elif(merger_choice == "njmerge"):
                    first_tree_dendropy = dendropy.Tree.get(path=first_tree, schema="newick")
                    second_tree_dendropy = dendropy.Tree.get(path=second_tree, schema="newick")
                    _,pair_tree_dendropy = njmergepair.run(node_distance_matrix, node_distance_matrix_taxlist, first_tree_dendropy, second_tree_dendropy)
                    # pair_tree_dendropy.is_rooted = False
                    # pair_tree_dendropy.collapse_basal_bifurcation(set_as_unrooted_tree=True)

            pair_tree_dendropy.resolve_polytomies()
            pair_tree_dendropy.update_bipartitions()
            pair_tree_dendropy.write(path=pair_tree_filename, schema="newick", suppress_rooting=True)

            subprocess.call(["raxmlng", "--evaluate", "--msa", pair_fasta, "--threads", "2", "--model", "GTR+G", "--tree", pair_tree_filename, "--prefix", pair_tree_filename])
            subprocess.call(["mv", pair_tree_filename + ".raxml.bestTree", pair_tree_filename])

if __name__ == "__main__":
    run_merger_and_raxml()
