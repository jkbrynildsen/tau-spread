# tau-spread

Code to reproduce all network analyses in Lubben et al. "LRRK2 kinase inhibition reverses G2019S mutation-dependent effects on tau pathology spread". Scripts are modified versions of those found in https://github.com/ejcorn/tau-spread.

## Requirements:
  - Tested with R 4.2.1
    - Requisite packages are listed in code/packages.R

## Tutorial

For a general overview of linear diffusion modeling in neurodegenerative disease, see tutorial in `example` folder here: https://github.com/ejcorn/tau-spread.

## Directory structure

Master branch contains 2 major folders:
  - `code` contains folders with scripts organized by analysis, i.e. `code/diffmodel` contains code that uses linear diffusion models to predict spread of protein through structural connectome.
  - `data` contains input structural data and pathology data.

## Input specification

The file `pipeline.R` is located in the main directory. This file will coordinate the sequential execution of all scripts within the code/ folder, generating all the figures in the paper and more. Custom specification of the following inputs at the top of `pipeline.R` is required:
  - `basedir`:  path to the main directory containing the `code` and `data` folders 
  - `matlab.path`: path to MATLAB binary. only needed for one supplemental analysis.
  - `injection.site`: vector of character names of brain regions in ABA to inject pathology into
  - `grps`: character vector containing the name of genotypes in data file to test. For our dataset, these were 'NTG' and 'G20'.
  - `treatments`: character vector containing the names of treatment groups in the data file to test. For our dataset, these were '0 MLi-2', '75 MLi-2', and '450 MLi-2'.
  - `opdir`: here you can add some extra custom label for your output folder given a particular configuration at the top of the script

## Questions, suggestions, comments?

Please contact Kate Brynildsen (jbryn ~`AT`~ seas.upenn.edu) with any questions regarding network analysis and code, and contact Mike Henderson (Michael ~`DOT`~ Henderson ~`AT`~ vai.org) with any questions regarding experiments and data.
