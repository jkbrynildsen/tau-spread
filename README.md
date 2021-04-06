# tau-spread

Code to reproduce all analysis in Henderson et al. 2021 ("Tau pathology spreads between anatomically-connected regions of the brain and is modulated by a LRRK2 mutation")

## Requirements:
  - Tested with R 3.6.1 
    - Requisite packages are listed in code/packages.R
  - MATLAB, tested with R2019a, in order to make network null models in the supplement

## Tutorial

For a general overview of linear diffusion modeling in neurodegenerative disease, see tutorial in `example` folder.

## Directory structure

Master branch contains 2 major folders:
  - `code` contains folders with scripts organized by analysis, i.e. `code/diffmodel` contains code that uses linear diffusion models to predict spread of protein through structural connectome, `code/G20vsNTG` deals with comparisons between G20 mice and NTG mice.
  - `data` contains input pathology data, look-up table to convert ABA and CNDR spaces

## Input specification

The file `pipeline.R` is located in the main directory. This file will coordinate the sequential execution of all scripts within the code/ folder, generating all the figures in the paper and more. Custom specification of the following inputs at the top of `pipeline.R` is required:
  - `basedir`:  path to the main directory containing the `code` and `data` folders 
  - `matlab.path`: path to MATLAB binary. only needed for one supplemental analysis.
  - `injection.site`: vector of character names of brain regions in ABA to inject pathology into
  - `injection.site.CNDR`: named vector that will be used to convert ABA to CNDR regions for injection sites only
  - `grps`: character vector containing the name of groups in data file to test. For our data set, these were 'NTG' and 'G20'.
  - `opdir`: here you can add some extra custom label for your output folder given a particular configuration at the top of the script

## Using individual scripts

The pipeline, when sequentially executed, will produce all analyses in the paper. However, some analyses are very time consuming (such as the `*traintest*.R` scripts) and it is desirable to be able to run individual scripts in isolation. In order to do so, execute `pipeline.R` down to line 25. Then you can open almost any script and execute it in isolation.

## Questions, suggestions, comments?

Please contact Eli Cornblath (Eli ~`DOT`~ Cornblath ~`AT`~ pennmedicine.upenn.edu) with any questions regarding network analysis and code, and contact Mike Henderson (hendm ~`AT`~ pennmedicine.upenn.edu) with any questions regarding experiments and data.
