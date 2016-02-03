PPI simulation work of Kacharoo et al., Science 348:921–925, 2015
========================
This repository contains the code and data that were used to carry out the protein-protein interaction simulations of this paper:
> A. H. Kachroo, J. M. Laurent, C. M. Yellman, A. G. Meyer, C. O. Wilke, E. M. Marcotte (2015). Systematic humanization of yeast genes reveals conserved functions and genetic modularity. Science 348:921–925.

It is provided as is, as a record of the research performed. The code in this repository is not actively maintained.

## Description of the analysis pipeline

The main setup folder for running the evolutionary simulation is the PPI_Setup directory.  To start a simulation, the whole directory can be copied as it stands.  Once the directory is completely copied, a simulation can be started with the command `python main.py`.  At that point, it will fail and show what a successful command line input will look like.

The code depends on many things.  First, biopython must be installed and importable with the `import Bio` command.  Next the location of Foldx is currently hardcoded as `/home/agm854/local/bin/foldx64Linux`.  Of course, that will need to be changed for your location.  The location of Foldx may be a command line option in the future.  

After completion of the simulations, the `analyze_ppi.py` script should be run to build the ancestral reconstructions and save the information gleaned from it.  To run the analysis, you will run `python analyze_ppi.py kept_mutants.txt start_structure.pdb`.  In addition to Biopython, the analysis requires numpy.

To get the full data files, you then need to run `interface_identity.py` with the command `python interface_identity.py kept_mutants.txt distance_cutoff start_structure.pdb ancestral_comparisons.txt`.  In addition to standard libraries, the script `interface_identity.py` requires numpy, prody (from the Bahar lab), and Biopython.  The final data file will be called `interface_identity.txt`.

The plots are all made with R.  All of the R scripts are in the PPI_analysis folder.  You will need to change the paths and have at least the survival, ggplot2, and latticeExtra libraries installed.  All of the data that we used in the manuscript is in the `data` folder and all of the scripts create figures in the `figures` folder.

