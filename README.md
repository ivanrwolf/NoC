# Node Classification - NoC
A node classification tool for gene regulatory networks (GRNs).

**Requirements:**
 * Linux, or unix like OS with support to conda package mananger
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

In order to provide a easy setup please use conda or miniconda.
Install instruction for conda https://docs.anaconda.com/anaconda/install/index.html.

With conda installed, you can install all dependecies with:
 >conda env create -f conda_NoC_environment.yml

Then, activate conda enviropnment:
 >conda activate NoC

**Usage:**

With conda environment activated you can simply run:
 >python grn.py [input]

**Positional arguments:**
 * input
 >Input file, currently only supports "ncol" format.

**Example:**
 >python grn.py 01_input_Scerevisiae_net.ncol

In order to exit from environment:
 >conda deactivate

If you use NoC please cite: [placeholder].
