# tau-spread

Code to reproduce all analysis in Henderson et al. 2021 ("Tau pathology spreads between anatomically-connected regions of the brain and is modulated by a LRRK2 mutation")

## Requirements:
  - Tested with R 3.6.1 
    - Requisite packages are listed in code/packages.R
  - MATLAB, tested with R2019a, in order to make network null models in the supplement

## Tutorial

For a general

## Directory structure

Master branch contains 2 major folders:
  - code/ contains folders, each of which contain scripts specific to certain analyses, i.e. ‘code/diffmodel’ contains code that uses linear diffusion models to predict spread of protein through structural connectome, 'code/G20vsNTG' deals with comparisons between G20 mice and NTG mice.

Using gene expression and linear dynamics of spread along the connectome, we attempt to predict the spatial distribution of experimentally observed pathology at 1, 3, and 6 months post injection.

## Input specification

The file ‘pipeline.R’ is located in the main directory. This file will coordinate the sequential execution of all scripts within the code/ folder, generating all the figures in the paper and more. Custom specification of the following inputs at the top of ‘pipeline.R’ is required:
  - basedir:  path to the main directory containing the 'code' and 'Data83018' folders 
  - matlab.path: path to MATLAB binary
  - opdir: name of output directory that contains all results, which will be housed in basedir
  - grps: character vector containing the name of groups in data file to test. For our data set, these were 'NTG' and 'G20'.

## Using individual scripts

## Questions, suggestions, comments?

Please contact Eli Cornblath (Eli ~`DOT`~ Cornblath ~`AT`~ pennmedicine.upenn.edu) with any questions regarding network analysis and code, and contact Mike Henderson (hendm ~`AT`~ pennmedicine.upenn.edu) with any questions regarding experiments and data.
