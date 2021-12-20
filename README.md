# Node Classification - NoC
A node classification tool for gene regulatory networks (GRNs).

**Requirements:**
 * Linux, or unix like OS with support to conda package mananger and bash shell
 * Weka 3.8.5
 * Python 3
 * Pandas
 * Numpy
 * scipy
 * igraph

**Installation:**

Download NoC files with:
 >git clone https://github.com/ivanrwolf/NoC
 
Change directory to the NoC folder:
 >cd ./NoC

Download Weka version 3.8.5 (https://waikato.github.io/weka-wiki/downloading_weka/) and paste it in NoC folder.
Alternatively you can edit NoC.sh and add the path to the "weka.jar" file manually.


In order to provide a easy setup please use conda or miniconda.
Install instruction for conda https://docs.anaconda.com/anaconda/install/index.html.

With conda installed, you can install all dependecies with:
 >conda env create -f conda_NoC_environment.yml

Then, activate conda enviropnment:
 >conda activate NoC

**Usage:**

With conda environment activated you can simply run:
 >bash NoC.sh [input]

**Positional arguments:**
 * input
 >Input file, currently only supports "ncol" format.

**Example:**
 >bash NoC.sh input_Scerevisiae_net.ncol

In order to exit from conda environment:
 >conda deactivate

If you use NoC please cite: Wolf, I.R., Simões, R.P. & Valente, G.T. Three topological features of regulatory networks control life-essential and specialized subsystems. Sci Rep 11, 24209 (2021). https://doi.org/10.1038/s41598-021-03625-w.
