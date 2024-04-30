#####################################################################
#
#  Random Forest integrative analysis tool
#
#                                    2019.01.30 Masahiro Ryo
#####################################################################

# To define the two functions for permutation-based random forest
#   1. myvarimp
#   2. RF_Hapfelmeier
#
#   INPUT: none
#   OUTPUT: none (This script is called by the main script: "SML_application.R")
#

  # install devtools
  tmp.install = require("devtools", character.only = TRUE)==FALSE
  if(tmp.install) install.packages("devtools")
  # install patchwork
  tmp.install = require("patchwork", character.only = TRUE)==FALSE
  if(tmp.install) devtools::install_github("thomasp85/patchwork")
  # install and require them
  package.list = c("ggplot2", "mlr", "patchwork", "party", "caret","foreach","doSNOW", "readxl", "mmpf")
  tmp.install = which(lapply(package.list, require, character.only = TRUE)==FALSE)
  if(length(tmp.install)>0) install.packages(package.list[tmp.install])
  suppressMessages(lapply(package.list, require, character.only = TRUE))

###-------------------------------------------------------------------------
### 1. myvarimp
###    function for calculating variable importance measure (modified by Masahiro Ryo)
###-------------------------------------------------------------------------
  myvarimp = function(object, mincriterion = 0, conditional = FALSE, 
                      pre1.0_0 = conditional, varID) {
    response = object@responses
    input = object@data@get("input")
    xnames = colnames(input)
    inp = initVariableFrame(input, trafo = NULL)
    y = object@responses@variables[[1]]
    if (length(response@variables) != 1) stop("cannot compute variable importance measure")
    CLASS = all(response@is_nominal)
    ORDERED = all(response@is_ordinal)
    if (CLASS) {
      error = function(x, oob) mean((levels(y)[sapply(x, which.max)] != y)[oob])
    }else {
      if (ORDERED) {
        error = function(x, oob) mean((sapply(x, which.max) != y)[oob])
      } else {
        error = function(x, oob) mean((unlist(x) - y)[oob]^2)}
    }
    perror = rep(0, length(object@ensemble))
    for (b in 1:length(object@ensemble)) {
      tree = object@ensemble[[b]]
      w = object@weights[[b]]
      w[w == 0] = 1
      oob = object@weights[[b]] == 0
      p = .Call("R_predict", tree, inp, mincriterion, -1L, PACKAGE = "party")
      eoob = error(p, oob)
      j = varID
      p = .Call("R_predict", tree, inp, mincriterion, as.integer(j), PACKAGE = "party")
      perror[(b - 1)] = (error(p,oob) - eoob)
    }
    return(MeanDecreaseAccuracy = mean(perror))
  }
  environment(myvarimp) = environment(varimp)
  
###-------------------------------------------------------------------------
### 2. Permutation-based random forest function (modified by Masahiro Ryo)
###   
###-------------------------------------------------------------------------
  RF_permutation = function(formula, data, nperm = 100, ntree = 100, ncore=2, varlim=-1, alpha=0.05)
    # formula: object of class "formula".
    # data: data frame containing the variables.
    # nperm: Number of permutation steps used for the permutation test.
    # ntree: Number of trees in the Random Forest.
    # ncore: Number of cores used for parallel computing
  {
    x.names = all.vars(formula)[-1]
    y.names = all.vars(formula)[1]
    terms. = terms.formula(formula)
    x.formula = attr(terms., "term.labels")
    y.formula = as.character(attr(terms., "variables"))[2]
    mtry = ceiling(sqrt(length(x.formula)))
    dat = subset(data, select = c(y.names, x.names))
    forest = party::cforest(formula, data = dat, controls = cforest_unbiased(mtry = mtry, ntree = ntree))
    obs.varimp = varimp(forest)
    obs.varimp[which(obs.varimp<0)] = 0
    perm.mat = matrix(10^5, ncol = length(x.names), nrow = nperm, dimnames = list(1:nperm, x.names))
    
    if(varlim == -1) varlim = length(x.formula)
    if(varlim > length(x.formula)) varlim = length(x.formula)
    x.names.selected = names(sort(obs.varimp, decreasing = T)[1:varlim])
    
    cl=makeCluster(ncore) #change to your number of CPU cores
    registerDoSNOW(cl)
    
    for (j in x.names.selected) {
      cat("\r", "Processing variable ", which(j == x.names.selected), " of ", length(x.names.selected)); flush.console()
      perm.dat = dat
      perm.mat[, j] = unlist(foreach(i = 1:nperm, .packages='party',.export="myvarimp") %dopar% {
        perm.dat[, j] = sample(perm.dat[, j]);
        myvarimp(party::cforest(formula, data = perm.dat, controls = cforest_unbiased(mtry = mtry, ntree = ntree)), varID = which(x.names == j))
      })
    }
    stopCluster(cl)
    
    p.vals = sapply(x.names, function(x) sum(perm.mat[, x] >= obs.varimp[which(x == x.names)]) / nperm)
    p.vals.bonf = p.vals * length(p.vals)
    
    if (any(p.vals < alpha)) {
      selection = names(p.vals)[which(p.vals < alpha)]
      mtry = ceiling(sqrt(length(selection)))
      forest = cforest(as.formula(paste(y.formula, ".", sep = " ~ ")), data = subset(dat, select = c(y.names, selection)),
                       controls = cforest_unbiased(mtry = mtry, ntree = ntree))
      p = p.vals[which(p.vals < alpha)]
      accuracy.fitting    = caret::postResample(predict(forest), subset(dat, select = y.names))
      accuracy.validation = caret::postResample(predict(forest, OOB=T), subset(dat, select = y.names))
      varimp.R2 = accuracy.fitting[2]*obs.varimp/sum(obs.varimp)
    }
    if (any(p.vals.bonf < alpha)) {             
      selection.bonf = names(p.vals.bonf)[which(p.vals.bonf < alpha)]                         
      mtry = ceiling(sqrt(length(selection.bonf)))
      forest.bonf = party::cforest(as.formula(paste(y.formula, ".", sep = " ~ ")), data = subset(dat, select = c(y.names, selection.bonf)),
                                   controls = cforest_unbiased(mtry = mtry, ntree = ntree))
      p.bonf = p.vals.bonf[which(p.vals.bonf < alpha)]
      varimp.bonf = varimp(forest.bonf)
      varimp.bonf[which(varimp.bonf<0)] = 0
      accuracy.fitting.bonf = caret::postResample(predict(forest.bonf), subset(dat, select = y.names))
      accuracy.validation.bonf = caret::postResample(predict(forest.bonf, OOB=T), subset(dat, select = y.names))
      varimp.R2.bonf = accuracy.fitting[2]*varimp.bonf/sum(varimp.bonf)
    }
    if (!any(p.vals < alpha)) {
      selection = c(); forest = c(); p = c();varimp = c();
      accuracy.fitting = c(); accuracy.validation = c(); varimp.R2 = c(); residual = c()
    }
    if (!any(p.vals.bonf < alpha)) {
      selection.bonf = c(); forest.bonf = c(); p.bonf = c(); varimp.bonf = c();
      accuracy.fitting.bonf = c(); accuracy.validation.bonf = c(); varimp.R2.bonf = c(); residual.bonf = c()
    }
#    Y = as.numeric(as.character(dat[,y.names]))
    # oob.error = ifelse(length(selection) != 0, mean((Y - as.numeric(as.character(predict(forest, OOB = T))))^2), mean((Y - ifelse(all(Y %in% 0:1), round(mean(Y)), mean(Y)))^2))
    # oob.error.bonf = ifelse(length(selection.bonf) != 0, mean((Y - as.numeric(as.character(predict(forest.bonf, OOB = T))))^2), mean((Y - ifelse(all(Y %in% 0:1), round(mean(Y)), mean(Y)))^2))
    cat("\n", "\r"); flush.console()
    return(list("forest" = forest, "forest.bonf" = forest.bonf, 
                "fitting" = accuracy.fitting, "fitting.bonf" = accuracy.fitting.bonf, 
                "validation" = accuracy.validation, "validation.bonf" = accuracy.validation.bonf,
                "p.values" = p.vals, "p.values.bonf" = p.vals.bonf,
                "varimp" =obs.varimp, "varimp.bonf"=varimp.bonf,  
                "varimp.R2"=varimp.R2, "varimp.R2.bonf"=varimp.R2.bonf))
  }
  
###-------------------------------------------------------------------------
### 3. is there any character variable?
###-------------------------------------------------------------------------
  
  check.character = function(data){
    data = data.frame(data)
    if(any(sapply(data,is.character))){
      id = which(sapply(data,is.character)==TRUE)
      warning("check the following variable(s) type: character?")
      warning(colnames(data)[id])
    }else{
    print("no character type detected")
    }
  }