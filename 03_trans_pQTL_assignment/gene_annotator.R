#############################################
## function to build gene annotator based on
## a given true positive set and features

gen.annotator.rf <- function(dat, pred, trc){
  
  ## 'dat'  -- data set containing all relevant annotations
  ## 'pred' -- set of features to be used for prediction
  ## 'trc'  -- column for true positive findings
  
  ## --> create split containing pseudo true positive examples <-- ##
  
  ## get all R2.counts with at least one true positive finding 
  ## (careful: may reassign one, if the actual lead SNP is not a VEP top consequence, but another SNP in the credible set)
  train.dat <- dat[ eval(as.name(trc)) == 1 ]$testing.flag
  ## keep only one 'pheno' per protein
  train.dat <- dat[ testing.flag %in% train.dat ]
  train.dat <- train.dat[ order(testing.flag, hgnc_symbol, pheno)]
  train.dat[, ind := 1:.N, by=c("testing.flag", "hgnc_symbol")]
  ## retain only unique examples
  train.dat <- train.dat[ ind == 1]
  
  ## print numbers
  cat("\n", nrow(train.dat), "observations for training\n")
  ## how many pseudo true positive examples
  cat("\n", table(train.dat[, trc, with=F])[2], "positive examples\n")
  
  #-------------------------------------------#
  ## --> build model using random forest <-- ##
  #-------------------------------------------#
  
  ## --> create test/validation sets <-- ##
  
  ## ensure that r2.groups are separated (including secondary signals) 
  ## protein not possible due to pleiotropic effects of some variants
  pqtl.sample        <- unique(train.dat[ eval(as.name(trc)) == 1, .(R2.group, hgnc_symbol)])
  cat("\n", nrow(pqtl.sample), "r2 group examples examples\n")
  
  ## do some renaming
  names(train.dat)   <- gsub("-", "_", names(train.dat))
  
  ## convert to data frame to ease application
  train.dat          <- as.data.frame(train.dat)
  
  ## allow parallel
  registerDoMC(20)
  
  ## build ten models, with different splits each
  pqtl.classifier.10 <- lapply(1:10, function(x){
    
    ## draw train and validation set
    set.seed(x)
    ii               <- sample(1:nrow(pqtl.sample), nrow(pqtl.sample)*.7)
    ii.train         <- pqtl.sample[ ii,]
    ii.test          <- pqtl.sample[-ii,]
    
    
    ## train the classifier
    pqtl.classifier  <- train(as.matrix(train.dat[ train.dat$R2.group %in% ii.train$R2.group, pred]), 
                              as.factor(train.dat[ train.dat$R2.group %in% ii.train$R2.group, trc]), 
                              method="rf", metric = "Kappa",
                              tuneGrid=expand.grid(mtry=2:6),
                              trControl= trainControl(method = "repeatedcv", number = 3, repeats = 2, 
                                                      sampling = "rose",
                                                      allowParallel = T))
    
    ## obtain prediction metrics in independent test set
    pqtl.pred.test   <- predict(pqtl.classifier, type = "prob", 
                                newdata = as.matrix(train.dat[ train.dat$R2.group %in% ii.test$R2.group, pred]))
    
    ## store confusion matrix with prob .5 as theshold
    conf.mat         <- confusionMatrix(as.factor(ifelse(pqtl.pred.test[, 2] > .5, 1, 0)), 
                                        as.factor(train.dat[ train.dat$R2.group %in% ii.test$R2.group, trc]))
    ## store roc curve
    pqtl.roc         <- roc(train.dat[ train.dat$R2.group %in% ii.test$R2.group, trc], 
                            pqtl.pred.test[, 2])
    
    ## return output
    return(list(pqtl.rf=pqtl.classifier, pqtl.pred=pqtl.pred.test, conf.mat=conf.mat, pqtl.roc=pqtl.roc))
    
  })
  
  ## return classifier
  return(pqtl.classifier.10)
}