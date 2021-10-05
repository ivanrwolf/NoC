# NoC
A node classification tool for gene regulatory networks (GRNs)

**Requirements:**
In order to provide a easy installation please use anaconda (https://docs.anaconda.com/anaconda/install/index.html) then run:
>conda env create -f conda_NoC_environment.yml

Then:
>conda activate NoC

In order to exit from environment:
>conda deactivate

**usage:**
>python grn.py [input]

**positional arguments:**
 * input
 >Input File, currently only supports "ncol" format.

**example:**
>python grn.py 01_input_Scerevisiae_net.ncol
