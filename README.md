# Bayesian-DNA-Barcode-Gap-Coalescent

A Bayesian implementation coded in Stan to a nonparametric Maximum Liklihood approach introduced here: https://github.com/jphill01/DNA-Barcode-Gap-Coalescent

## How to Run

``Analysis.R`` contains one function, `run_DNA_barcode_gap_analysis_by_marker()`, which prompts users interactively to select the folder containing marker subfolders (*e.g.*, COI-3P, COI-5P, CYTB) and the `DN_barcode_gap.stan` program located on their desktop. 

Data needed to run ``Analysis.R`` can be found in the ``Markers`` folder.
