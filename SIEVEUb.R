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

family_bootstrap_svm <- function(data, factor, family, training_division=0.5, niter=10, 
                                 downsample_families=NA, oldmean=NA) {
	# takes data as a set of features (rows) and examples (columns) such as
	#       from microrarry data
	# family is a vector of values for family calls that corresponds to the items in factor. This identifies
  #       groups of examples that are similar to each other- for example, when examined using traditional
  #       methods.
	# A corresponding set of factors must be included which are named as
	#       the columns of the data and are binary for positive and negative examples in the dataset
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

	
	# given that the selection of training families is random, we should run this process a number
	#    of times to ensure that we're getting a good sampling of these divisions.
	#    Thus we'll run niter times
	for (i in 1:niter) {
		# downsample families. This provides an optional upper threshold on the size
	  #     of individual families. This is mainly to test whether the presence of large 
	  #     families could skew results (it doesn't seem to).
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
		# to do this we select entire families for training division instead of individual examples.
		#       This process ensures that 'trivial' relationships in family members don't artficially
		#       inflate the performance. That would happen by selection of one family member for the
		#       training set, and another example from the same family going in to the testing set.
		#       family-wise cross-validation prevents this from happening and provides a conservative
		#       estimate of performance.
		trx = sample(families, length(families)*training_division)
		tex = families[which(!families %in% trx)]
		
		# now use the families to select the corresponding examples
		training = names(these_families)[which(these_families %in% trx)]
		testing = names(these_families)[which(these_families %in% tex)]
		
		#browser()
		
		# train a model based on the training set and...
		thismodel = try(svm(x=t(data[,training]), factor[training], probability=T), TRUE)
		if (class(thismodel)=="try-error") next
		
		# evaluate it's performance on the held out testing set
		thesepredictions = try(predict(thismodel, t(data[,testing]), probability=T), TRUE)
		if (class(thesepredictions)=="try-error") next
		
		models = c(models, list(thismodel))
			
		# calculate an ROC AUC from the prediction as the performance
		class_auc = try(auc(factor[testing], attr(thesepredictions, "probabilities")[,1]), TRUE)
		if (class(class_auc) == "try-error") class_auc = NA
		
		prob[testing,i] = attr(thesepredictions, "probabilities")[testing,1]
		
		cat(".")
		
		res = c(res, class_auc)
	}
	
	# There was an issue with this process. It seems that when the model loses it's ability
	#  to disciriminate between positive and negative examples the distribution of proteins in
	#  families starts producing highly odd results skewed in favor of some families
	#  and thus inflating the ability of the model to predict.
	# This issue has been fixed by normalizing the probabilities from each iteration prior to averaging.
	
	#if (oldmean) prob[,(i+1)] = sapply(rownames(prob), function (r) max(prob[r,1:niter], na.rm=T))
	#browser()
	if (is.na(oldmean)) {
	  prob[,1:niter] = sapply(1:niter, function (i) prob[,i]/max(prob[,i], na.rm=T))
	}
	prob[,(i+1)] = sapply(rownames(prob), function (r) mean(prob[r,1:niter], na.rm=T))
	
	#browser()
	
	overall = auc(factor, prob[,"mean"])
	
	#browser()
	return(list(auc=overall, aucs=res, probs=prob, models=models))
}

svmrfeFeatureRanking = function(x,y, stepamount=0.5) {
  # where x is the data - columns are features
  #       y is a class factor for the examples
  n = ncol(x)
  survivingFeaturesIndexes = seq(1:n) 
  featureRankedList = colnames(x)
  scoresRankedList = vector(length=n) 
  rankedFeatureIndex = n
  step_length = length(featureRankedList)
  
  results = list()
  i = 0 
  while(step_length>2) {
    i = i + 1
    
    #train the support vector machine
    svmModel = svm(x[, featureRankedList], y, scale=F, probability=T)
    #return(svmModel)
    #compute the weight vector
    w = t(svmModel$coefs)%*%svmModel$SV
    #w = svmModel$coefs
    #compute ranking criteria
    rankingCriteria = w * w
    #rank the features
    ranking = sort(rankingCriteria, index.return = TRUE)
    
    stepx = as.integer(length(ranking$ix)*stepamount)
    
    # remove the lowest portion of the ranked features
    removal = ranking$ix[1:stepx]
    scores = w[ranking$ix][1:stepx]
    
    #thesepredictions = predict(svmModel, x[,survivingFeaturesIndexes], probability=T)
    
    #step_auc = try(auc(y, attr(thesepredictions, "probabilities")))
    step_auc = auc(y, unlist(svmModel$decision.values[,1]))
    
    step_length = length(featureRankedList)
    
    scoresRankedList = scores
    featureRankedList = colnames(rankingCriteria)[ranking$ix[(stepx+1):length(ranking$ix)]]
    
    #update weights
    #scoresRankedList[rankedFeatureIndex:(rankedFeatureIndex-length(scores))] = scores
    
    #update feature ranked list
    #featureRankedList[rankedFeatureIndex:(rankedFeatureIndex-length(removal))] = survivingFeaturesIndexes[removal]
    
    #eliminate the features with smallest ranking criterion 
    #survivingFeaturesIndexes = survivingFeaturesIndexes[-removal]
    #rankedFeatureIndex = rankedFeatureIndex-length(removal)
    
    cat("Remaining features: ", length(featureRankedList), "\n")
    #browser()
    
    results[[i]] = list(nfeatures=step_length, auc=step_auc, features=featureRankedList, weights=scoresRankedList)
  }
  return (results)
}

svmrfeCrossFeatureRanking = function(x, y, family, stepamount=0.5, training_division=0.5, 
                                     niter=10, fbs_niter=10, downsample_families=NA) {
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
                                    training_division=training_division, niter=fbs_niter,
                                    downsample_families = downsample_families)
    
    # now calculate weights for features
    for (iter in 1:niter) {
      # randomly sample the feature space at the stepamount
      these_features = sample(survivingFeaturesIndexes, length(survivingFeaturesIndexes)*stepamount)
      #train the support vector machine
      svm_cross = try(family_bootstrap_svm(t(x[,these_features]), y, family, training_division=training_division, 
                                           niter=fbs_niter, downsample_families = downsample_families))
      if (class(svm_cross)=="try-error") next
      
      # now assign
      weights[iter,these_features] = rep(svm_cross$auc, length(these_features))
      cat(".")
    }
    
    # trim the weights to only those features still in consideration
    weights = weights[,survivingFeaturesIndexes]
    w = colMeans(weights, na.rm=T)
    
    step_auc = svm_base$auc
    # testing this
    #step_auc = mean(svm_base$aucs)
    
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
    
    results[[i]] = list(nfeatures=length(survivingFeaturesIndexes), auc=step_auc, aucs=svm_base$aucs, features=featureRankedList, weights=scoresRankedList)
    
    #eliminate the features with smallest ranking criterion 
    survivingFeaturesIndexes = survivingFeaturesIndexes[-removal]
    
  }
  return (results)
}

family_counts = function(data, families, these_fact) {
  # first we make sure we get rid of multiple counts from the same protein
  data = as.matrix(data)
  data[which(data>0)] = 1
  
  results = matrix(ncol=5, nrow=ncol(data))
  rownames(results) = colnames(data)

  for (i in 1:ncol(data)) {
    pcount = 0
    pscore = 0
    ncount = 0
    nscore = 0
    pos = 0
    neg = 0
    
    for (fam in unique(families)) {
      famil = names(which(families==fam))
      #print(fam)
      
      # we'll use a weighted score for each kmer in families
      count = sum(data[famil,i])
      score = count/length(famil)
      
      # this assumes that all in the family are labeled the same way
      # if they're not then that might be an issue overall
      if (these_fact[famil][length(famil)]=="positive") {
        pcount = pcount + count
        pscore = pscore + score
        pos = pos + 1
        #print(c(ncol(data), count, length(famil)))
        #browser()
      } 
      else {
        ncount = ncount + count
        nscore = nscore + score
        neg = neg + 1
      }
    }
    pscore = pscore/pos
    nscore = nscore/neg
    results[i,] = c(pcount, pscore, ncount, nscore, pscore-nscore)
  }
  return(results)
}

