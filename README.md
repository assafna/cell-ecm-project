# cell-ecm-project
This repository includes the code written for the analysis and quantification of long-range cell-cell mechnaical communication as presented in the paper:
https://www.biorxiv.org/content/10.1101/2020.07.30.223149v1


[figures.py](figures.py) includes all figures generated for the paper by order.

[fiji_scripts](fiji_scripts) includes all scripts used in Fiji software.

[fiber_density](fiber_density) folder includes all methods for analyses made based on quantification of the fibrin channel for both simulations and experiments.

[libs](libs) folder includes all methods used during the analysis process, including configurations and paths. Simulations-specific and experiments-specific methods are located in the inner folders.

[methods](methods) folder includes all methods used in the pre-processing such as creating new cell pairs, fake pairs, tracking of cells, printing data and computing the normalization parameters for z-score calculation.

[plotting](plotting) folder mainly include a modified version that uses plotly lib for creating and saving figures.

[data](data) folder includes the quantified windows of cell pairs for simulations and experiments in CSV format.