# SIEVE-Ub
Code and data for the SIEVE-Ub prediction algorithm for prediction of E3 ubiquitin ligase mimics

Accompanies the preprint: https://peerj.com/preprints/27292/

Files:
KmerFeatures.py : Python code that transforms a protein fasta file into kmer features suitable for use
                  with. the SIEVEUb.R code.
SIEVEUb.R : R code to train and validate models based on kmer features. Includes code for family-wise
	    cross-validation and family-wise scoring of overrepresented kmers.
SIEVEUbApply.R : an R script that allows reading in a set of features output from KmerFeatures and applying
            the predictive model from the paper to output predictions for each protein.
SIEVEUbCode.Rmd : RMarkdown file for doing the analysis described in the paper. This RMarkdown has code for
                  regenerating most of the figures, tables, and supplemental figures and tables from the paper.
make_prediction.csh : shell script that will provide predictions of ubiquitin ligase given an input fasta file.

data/1032069.3.PATRIC.faa : example PATRIC genome fasta file from a single genome
data/FamiliesConservative.txt : list of example proteins used for positive and negative examples, their annotations,
                    and identified sequence families
data/PATRIC.genomes : list of genomes from PATRIC used in the paper
data/PATRIC_secretion_systems.txt : Type III, IV, and VI secretion system components identified in each PATRIC genome used
data/SIEVEUb_best_model.Rdata : Rdata file containing the best Ub ligase model
data/SIEVEUb_minimal_model.Rdata : Rdata file containing the minimal Ub. ligase model with 10 features

This material was prepared as an account of work sponsored by an agency of the
United States Government.  Neither the United States Government nor the United
States Department of Energy, nor Battelle, nor any of their employees, nor any
jurisdiction or organization that has cooperated in the development of these
materials, makes any warranty, express or implied, or assumes any legal
liability or responsibility for the accuracy, completeness, or usefulness or
any information, apparatus, product, software, or process disclosed, or
represents that its use would not infringe privately owned rights.

Reference herein to any specific commercial product, process, or service by
trade name, trademark, manufacturer, or otherwise does not necessarily
constitute or imply its endorsement, recommendation, or favoring by the United
States Government or any agency thereof, or Battelle Memorial Institute. The
views and opinions of authors expressed herein do not necessarily state or
reflect those of the United States Government or any agency thereof.

		 PACIFIC NORTHWEST NATIONAL LABORATORY
			      operated by
				BATTELLE
				for the
		   UNITED STATES DEPARTMENT OF ENERGY
		    under Contract DE-AC05-76RL01830
