# cgMLST_Cornell_FSL_script

These scripts were developed for running a cgMLST analysis offline using cgMLST schemes. We developed this tool with the Listeria monocytogenes cgMLST scheme and rules in mind. 
More information regarding the L. monocytogenes can be found in the Moura et al., 2016 paper. 
To use the tool, one will need to edit all 3 files to match their system and will need to download the cgMLST database.


# at_matcher_per_allele.sh

Place this file in a directory named "lmo-mlst".
Create another directory and name it "cgMLST". Inside cgMLST, create two other directories and name them "dbs" and "processed". 
"dbs" is where you will place the cgMLST database scheme. "processed" is where your output files will appear after you run the tool.

Inside "at_matcher_per_allele.sh, change the line: work="/PATH_TO/cgMLST", so it matches the path to your newly created "cgMLST" directory.

This script requires that BLAST, MUSCLE and PARALLEL are installed in your system and accessible.

# cgMLST_loop.sh

Make sure you have all assemblies in the same directory. Then run:

sh cgMLST.sh <inpath>

Where <inpath> is the directory where your assemblies are located.
The outputs will be located in the cgMLS/processed/ directory.

# hier_clust_v2.R

This is an R script (requires R) that will look into your "results" files and compare each pair of files to identify the number of allelic differences between each pair.

You need to modify the path to match your system.
