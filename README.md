# HbA1c Cutoff Analysis

Author: Sarita Perez

License: MIT

This repository contains a SAS 9.4 analysis workflow for evaluating HbA1c as a screening test for type 2 diabetes. The analysis simulates an adult population calibrated to published national estimates and examines how different HbA1c cutoffs affect screening performance.

## Overview

The main script, HbA1c_Cutoff_Analysis.sas, performs the following:

- simulates adult population health data with age, BMI, race/ethnicity, diabetes status, HbA1c, and fasting plasma glucose
- evaluates screening performance at candidate HbA1c cutoffs of 5.7%, 6.0%, 6.3%, and 6.5%
- reports sensitivity, specificity, positive predictive value, negative predictive value, and screen-positive rate
- generates ROC curves and Youden's Index for cutoff selection
- fits logistic regression models with and without age/BMI adjustment
- examines disparities in diabetes prevalence and screening performance by race/ethnicity

## Repository Structure

- [README.md](README.md) — project description and usage notes
- [SAS Program/HbA1c_Cutoff_Analysis.sas](SAS%20Program/HbA1c_Cutoff_Analysis.sas) — main SAS analysis workflow
- [Tables/](Tables/) — folder for generated tables and report outputs

## Requirements

- SAS 9.4 or later
- SAS/STAT for logistic regression procedures
- SAS/GRAPH or ODS graphics support for plots and report output

## How to Run

1. Open HbA1c_Cutoff_Analysis.sas in SAS Studio, SAS Enterprise Guide, or a local SAS installation.
2. Update the output file paths in the program so they point to locations on your machine.
3. Run the full program.

## Outputs

The script produces:

- descriptive summaries and frequency tables
- ROC analysis and Youden-based cutoff comparison
- logistic regression results
- PDF and RTF output reports
- RTF-format table files stored in [Tables/](Tables/) for sharing and documentation

## Notes

- This project uses simulated data rather than real patient-level records.
- The results are intended for methodological demonstration, screening analysis, and educational use.

