# TreeMerge
---

This repository contains scripts for generating pair-wise merged tree files and a new entry point to [TreeMerge](https://github.com/ekmolloy/treemerge) with additional input files.

Here is an overview of each script. For the full command line reference, you can also run the python files with the --help flag.

### get\_merge\_list.py
This file is used to generate a minimum spanning tree as well as to generate a file which defines which pairs of trees must be merged together in order for TreeMerge to ingest them.
'''
python get\_merge\_list.py --starting-tree <input tree> --output-prefix <output prefix> -- <space separated list of constraint trees>
'''

### setup\_merger.py
This file is used to compute the induced starting trees needed for the pair-wise mergers
'''
python setup\_merger.py --starting-tree <input tree> --files-needed <output named files_needed from get_merge_list.py> --output-prefix <output prefix> --merger-choice njmerge --guide-choide induced
'''

### run\_merger.py
This file is used to run the pair-wise mergers and create pair-wise merged trees
'''
python run\_merger.py --starting-tree <input tree> --files-needed <output named files_needed from get_merge_list.py> --output-prefix <output prefix> --merger-choice njmerge --guide-choice induced
'''

### treemerge.py
This file is the original treemerge that was modified to take in a predetermined minimum spanning tree. It expects the pair-wise trees to already have been computed on the input minimum spanning tree with ther branch lengths already estimated.
'''
python treemerge.py -s <input tree> -m <input matrix> -x <input matrix taxon list> -o <output prefix> -p <paup binary> -w <working directory> --mst <minimum spannig tree> -t <space separated list of constrait trees>
'''
