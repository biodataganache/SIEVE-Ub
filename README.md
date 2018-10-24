Code to accompany the paper (currently preprint): https://peerj.com/preprints/27292/

# SIEVE-Ub
Code and data for the SIEVE-Ub prediction algorithm for prediction of E3 ubiquitin ligase mimics

1. KmerFeatures.py: takes a protein fasta file and produces a feature file based on parameters. Required biopython.
2. SIEVEUb.R: R functions for SVM application, cross-validation, and feature selection. Requires e1071 and pROC libraries.
3. SIEVEUbCode.Rmd: R markdown file that shows how the code is used with the data to evaluate the approach, generate
                    the final model and apply it to a new genome to make predictions.
4. Util/: Python functions required by KmerFeatures.py
5. data/: relevant data files for examples
