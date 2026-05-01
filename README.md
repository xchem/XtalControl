# XtalControl v0.2
---
## Overview

XtalControl is a simple program made to analyse XChemExplorer (XCE) output files and return HTML reports on [tests of protein crystal survival of chemistry reagents](https://www.protocols.io/view/assessing-protein-crystal-tolerance-of-chemistry-r-14egn516yg5d/v2).  The program collates and averages repeats of the same compound or reagent and compares them to a set of editable benchmarks to return clear pass-fail results for each reagents feasability in [crude reaction mixture soaking experiments](https://www.nature.com/articles/s42004-020-00367-0) (i.e. does the reagent destroy the crystal or ruin the diffraction data). It provides a colour-coded and simple to understand assessment of crystal issues, mounting success, diffraction success, and diffraction data quality.

Successful diffraction of a crystal in a protein-ligand soak produces a clearly defined result that is usually easy to statistically analyse (e.g. *"we soaked 500 ligands, mounted and diffracted 400 crystals, and got 40 fragment hits"*). For other experiments, especially across multiple repeats, it is challenging and time-consuming to scroll through 150+ rows and 100+ columns to get an accurate feel for the experimental results. 

It would be fairly simple to repurpose the "SoakCode" logic in this program to assess other experiments (e.g. for different chemical series).

---
## Requirements

- [XChemExplorer](https://github.com/tkrojer/XChemExplorer)
    - A .csv export of your XCE database 
- A .csv platemap following the format shown in 'control-plate-map.csv'
- Python Environment containing:
    - [pandas](https://pandas.pydata.org/)
    - [numpy](https://numpy.org/)
    - [natsort](https://github.com/SethMMorton/natsort)
