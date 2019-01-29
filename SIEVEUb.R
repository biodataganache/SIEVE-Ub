# functions to perform training and evaluation of a model

library(e1071)
library(pROC)

bootstrap_svm <- function(data, factor, training_division=0.5, niter=10) {
	# takes data as a set of features (rows) and examples (columns) such as
	#       from microrarry data
	# A corresponding set of factors must be included which are named as
	#       the columns of the data and are binary
	# Training_division is how much is held out for training versus testing
	#       on each iteration
	# niter is the number of times to repeat the bootstrapping
	res = c()
	prob = matrix(nrow=length(factor), ncol=niter+1)
	rownames(prob) = names(factor)
	colnames(prob) = c(1:niter, "mean")
	
	for (i in 1:niter) {
		training = sample(names(factor), length(factor)*training_division)
		testing = names(factor)[which(!names(factor) %in% training)]
	
		thismodel = svm(x=t(data[,training]), factor[training], probability=T)
		thesepredictions = predict(thismodel, t(data[,testing]), probability=T)
				
		class_auc = auc(factor[testing], attr(thesepredictions, "probabilities")[,1])
		
		prob[testing,i] = attr(thesepredictions, "probabilities")[testing,1]
		
		res = c(res, class_auc)
	}
	prob[,(i+1)] = sapply(rownames(prob), function (r) mean(prob[r,1:niter], na.rm=T))
	overall = auc(factor, prob[,"mean"])
	
	return(list(auc=overall, aucs=res, probs=prob))
}

validate_svm <- function(train_data=NA, train_model=NA, test_data=NA, train_factor=NA, test_factor=NA, fill.missing.values=T) {
	# takes a training data set and an independent test data set
  # alternately takes an already trained model and applies it to a new feature set
	# constructs an SVM from the training data then applies it to the 	
	#       testing data and reports an AUC and the predictive model
	
	if (fill.missing.values) {
		test_data[which(is.na(test_data))] = 0
		train_data[which(is.na(train_data))] = 0
	}
	
  vids = rownames(test_data)
  thismodel = train_model
	if (is.na(train_model)) {
	  # get valid feature ids for the test set
	  vids = rownames(test_data)[rownames(test_data) %in% rownames(train_data)]
	  thismodel = svm(x=t(train_data[,names(train_factor)]), train_factor, probability=T)
	}
	
	thesepredictions = predict(thismodel, t(test_data[vids,]), probability=T)
	
  auc = 0
	if (!is.na(test_factor)) {
		# this should be a two factor list
		auc = auc(test_factor, attr(thesepredictions, "probabilities")[,1])
	}
	return(list(test_auc=auc, model=thismodel, predictions=attr(thesepredictions, "probabilities")))
}

family_bootstrap_svm <- function(data, factor, family, training_division=0.5, niter=10, downsample_families=NA) {
	# takes data as a set of features (rows) and examples (columns) such as
	#       from microrarry data
	# family is a vector of values for family calls that corresponds to the items in factor
	# A corresponding set of factors must be included which are named as
	#       the columns of the data and are binary
	# Training_division is how much is held out for training versus testing
	#       on each iteration
	# niter is the number of times to repeat the bootstrapping
  # downsample_families: if not NA should be an integer that limits family size to that number or smaller
  #                      of members. Larger families are randomly downsampled to that size.
	res = c()
	prob = matrix(nrow=length(factor), ncol=niter+1)
	rownames(prob) = names(factor)
	colnames(prob) = c(1:niter, "mean")
	models = list()
	
	cat("Running Family-wise bootstrap SVM ", niter, " times\n")
	
	families = unique(sort(family))

	for (i in 1:niter) {
		# downsample families
		these_families = family
		if (!is.na(downsample_families)) {
		  these_families = c()
		  
		  for (j in families) {
		    this_one = family[which(family == j)]
		    if (length(this_one)>downsample_families) this_one = sample(this_one, downsample_families)
		    these_families = c(these_families, this_one)
		  }
		}
	
		#ensure that examples are taken out as families
		trx = sample(families, length(families)*training_division)
		tex = families[which(!families %in% trx)]
		
		training = names(these_families)[which(these_families %in% trx)]
		testing = names(these_families)[which(these_families %in% tex)]
		
		
		thismodel = try(svm(x=t(data[,training]), factor[training], probability=T), TRUE)
		if (class(thismodel)=="try-error") next
		thesepredictions = try(predict(thismodel, t(data[,testing]), probability=T), TRUE)
		if (class(thesepredictions)=="try-error") next
		
		models = c(models, list(thismodel))
			
		class_auc = try(auc(factor[testing], attr(thesepredictions, "probabilities")[,1]), TRUE)
		if (class(class_auc) == "try-error") class_auc = NA
		
		prob[testing,i] = attr(thesepredictions, "probabilities")[testing,1]
		
		cat(".")
		
		res = c(res, class_auc)
	}
	prob[,(i+1)] = sapply(rownames(prob), function (r) mean(prob[r,1:niter], na.rm=T))
	overall = auc(factor, prob[,"mean"])
	
	return(list(auc=overall, aucs=res, probs=prob, models=models))
}

svmrfeCrossFeatureRanking = function(x, y, family, stepamount=0.5, training_division=0.5, niter=10, fbs_niter=10) {
  # where x is the data - columns are features
  #       y is a class factor for the examples
  n = ncol(x)
  survivingFeaturesIndexes = seq(1:n) 
  featureRankedList = vector(length=n)
  scoresRankedList = vector(length=n) 
  rankedFeatureIndex = n
  
  cat("Starting features: ", length(survivingFeaturesIndexes))
  
  results = list()
  i = 0 
  while(length(survivingFeaturesIndexes)>2) {
    i = i + 1
    
    # New: we need to produce the weights for the features here
    # For simplicity we pretend that all the features are here when we make
    #     the matrix. The only features that will get used are the survivingFeaturesIndexes
    #     and we trim the weights below.
    weights = matrix(NA, nrow=niter, ncol=ncol(x))
    
    # get the 'base' performance for this set of features
    svm_base = family_bootstrap_svm(t(x[,survivingFeaturesIndexes]), y, family, 
                                    training_division=training_division, niter=fbs_niter)
    
    # now calculate weights for features
    for (iter in 1:niter) {
      # randomly sample the feature space at the stepamount
      these_features = sample(survivingFeaturesIndexes, length(survivingFeaturesIndexes)*stepamount)
      #train the support vector machine
      svm_cross = try(family_bootstrap_svm(t(x[,these_features]), y, family, training_division=training_division, niter=fbs_niter))
      if (class(svm_cross)=="try-error") next
      
      # now assign
      weights[iter,these_features] = rep(svm_cross$auc, length(these_features))
      cat(".")
    }
    
    # trim the weights to only those features still in consideration
    weights = weights[,survivingFeaturesIndexes]
    w = colMeans(weights, na.rm=T)
    
    step_auc = svm_base$auc
    
    #compute ranking criteria
    rankingCriteria = w * w
    #rank the features
    # ranking works on the w list of scores
    ranking = sort(rankingCriteria, index.return = TRUE)
    dividing_point = as.integer(length(ranking$ix)*stepamount)
    
    # remove the lowest portion of the ranked features
    # ranking$ix is to the w list of scores
    removal = ranking$ix[1:dividing_point]
    
    # we are using "keeping" to grab the actual labels from the original list of features
    keeping = survivingFeaturesIndexes[ranking$ix[dividing_point:length(survivingFeaturesIndexes)]]
   
    cat("\nRemaining features: ", length(survivingFeaturesIndexes), " AUC: ", step_auc)
    
    
    # scores are from all remaining features
    # these scores are in order from smallest to largest
    scores = w[ranking$ix]
    
    #update weights
    scoresRankedList = scores
    
    #update feature ranked list
    featureRankedList = colnames(x)[keeping]
    
    #eliminate the features with smallest ranking criterion 
    survivingFeaturesIndexes = survivingFeaturesIndexes[-removal]
    
    results[[i]] = list(nfeatures=length(survivingFeaturesIndexes), auc=step_auc, features=featureRankedList, weights=scoresRankedList)
    
  }
  return (results)
}

svmCrossFeatureSensitivity = function(datamat, class_fact, families, training_division=0.5, feature_division=0.5, niter=1000) {
  n = ncol(datamat)
   
  # for this the objective is to capture information about the impact of features on prediction
  #     of examples and the impact of proteins in the training set on prediction of examples in
  #     the testing set. This information can tell us two main things: 1) for a given example
  #     what are the features most important to successfully predicting it, and 2) which other
  #     examples are the most responsible for its prediction- that is what are the most similar
  #     other examples?
  # This is a bit like RFE except without the elimination. We are only interested in doing it a
  #     bunch of times to build up good numbers on these relationships.
  
  # First we need to create results matrices
  
  protein_results = matrix(0, nrow=nrow(datamat), ncol=nrow(datamat), dimnames=list(rownames(datamat), rownames(datamat)))
  # give the counts a small value so we don't get divide by 0 NaNs
  protein_counts = matrix(0.01, nrow=nrow(datamat), ncol=nrow(datamat), dimnames=list(rownames(datamat), rownames(datamat)))
  feature_results = matrix(0, nrow=nrow(datamat), ncol=ncol(datamat), dimnames=dimnames(datamat))
  feature_counts = matrix(0.01, nrow=nrow(datamat), ncol=ncol(datamat), dimnames=dimnames(datamat))
  
  for (i in 1:niter) {
    # randomly select a set of features to test
    these_features = sample(colnames(datamat), ncol(datamat)*feature_division)
    
    # we need to 
    # 1) select a set of examples for training and then for testing
    # 2) train the svm on the training set
    # 3) predict on the testing set
    # 4) appropriately update the results matrices with the prediction probabilities
    
    # pulled from family_bootstrap_svm
    #ensure that examples are taken out as families
    trx = sample(1:max(families), max(families)*training_division)
   
    training = names(families)[which(families %in% trx)]
    testing = names(families)[which(!families %in% trx)]
    
    thismodel = try(svm(x=datamat[training,these_features], class_fact[training], probability=T), FALSE)
    if (class(thismodel)=="try-error") next
    thesepredictions = try(predict(thismodel, datamat[testing,these_features], probability=T), FALSE)
    if (class(thesepredictions)=="try-error") next
    
    predictions = attr(thesepredictions, "probabilities")[,1]
    
    #browser()
    
    protein_results[names(predictions), training] = protein_results[names(predictions), training] + predictions
    protein_counts[names(predictions),training] = protein_counts[names(predictions),training] + 1
    
    feature_results[names(predictions),these_features] = feature_results[names(predictions),these_features] + predictions
    feature_counts[names(predictions),these_features] = feature_counts[names(predictions),these_features] + 1
    
    #browser()

  }
  results = list(protein_results=protein_results/protein_counts, 
                 feature_results=feature_results/feature_counts, 
                 protein_counts=protein_counts, feature_counts=feature_counts)
  return (results)
}

# script to run SIEVEUb on example data
SIEVEUb = function() {
  # To generate features from a fasta file use the following command
  #    from the SIEVEServer code
  # SIEVEFeatures.py -f [fasta file] -o [output base] -m gist -s k:14:::reduced_alphabet_0
  
  
  # read in data file about examples (class and sequence family)
  ubex_classes = read.table("FamiliesConservative_10_04_2013.txt", sep="\t", header=1, row.names=1)
  ubex_fact = factor(x=ubex_classes[,3], labels=c("positive", "negative"))
  ubex_families = ubex_classes[,1]
  
  # this is a matrix of features from examples
  ubex_k14red0 = read.table("ubligase_k14red0.txt", sep="\t", row.names=1, header=1, stringsAsFactors=F)
  # filter out those features with minimal representation in examples
  ubex_k14red0r1 = ubex_k14red0[,which(sapply(colnames(ubex_k14red0), function (c) sum(ubex_k14red0[,c] != 0)>9))]
  
  # do a 100-fold cross validation
  ubex_k14red0r1_cv1 = family_bootstrap_svm(t(as.matrix(ubex_k14red0r1[names(ubex_fact),])), ubex_fact, ubex_families, niter=100)
  ubex_k14red0r1_cv1$auc
  
  # do a recursive feature elimination
  ubex_k14red0r1_rfe = svmrfeFeatureRanking(ubex_k14red0r1[ubex_fact,], ubex_fact, 0.25)
  
  #train model and apply it to salmonella genome
  salm_k14red0 = read.table("salmonella_lt2_k14red0.txt", sep="\t", row.names=1, header=1, stringsAsFactors=F)
  salm_k14red0r1 = salm_k14red0[,colnames(ubex_k14red0r1)]
  
  results = validate_svm(t(ubex_k14red0r1), t(salm_k14red0r1), ubex_fact)
  return(results)
} 


