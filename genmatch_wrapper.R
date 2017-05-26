library(foreign)
library(rgenoud)
library(Matching)
library(rbounds)

multiple.genmatch = function(dataset, outcome.label, match.formula, random.runs=1, pop.size=100, max.generations=50, wait.generations=20, thread.count=4) {
  outcome.vec = dataset[, outcome.label]
  var.labels = all.vars(match.formula)
  treatment.label = var.labels[1]
  treatment.vec = dataset[, treatment.label]
  covariates.labels = var.labels[2:length(var.labels)]
  
  covariates = c()
  for (i in c(1:length(covariates.labels))) {
    covariates = cbind(covariates, dataset[, covariates.labels[i]])
  }
  covariates = cbind(data.matrix(covariates))
  
  max.AMsmallest.p.value = -1
  current.best.result = list()
  for (i in c(1:random.runs)) {
    genmatch.result <- GenMatch(Tr=treatment.vec, X=covariates, BalanceMatrix=covariates, 
                                estimand="ATT", M=1, pop.size=pop.size, max.generations=max.generations, 
                                wait.generations=wait.generations, cluster=rep('localhost', thread.count))
    balance <- MatchBalance(
      match.formula,
      data=dataset,
      match.out=match.result,
      nboots=500
    )
    if (balance$AMsmallest.p.value > max.AMsmallest.p.value) {
      match.result <- Match(Y=outcome.vec, Tr=treatment.vec,
                            X=covariates, estimand="ATT", Weight.matrix=genmatch.result)
      current.best.result = list(genmatch.result, match.result, balance)
      max.AMsmallest.p.value = balance$AMsmallest.p.value
    }
  }
  return(current.best.result)
}

# use example:

result = multiple.genmatch(
  dataset,
  "outcome.label",
  treatment ~ covariate1 + covariate2 + covariate3,
  random.runs = 10,
  pop.size = 150,
  max.generations = 50,
  wait.generations = 25,
  thread.count = 8
)

# thread.count should usually be the logical CPU count (on a Mac, usually double the physical CPU core count:
# for a dual-core CPU, thread.count = 4; for a quad-core, thread.count = 8.)

# All arguments including random.runs and after are omittable.
# If omitted, they will take on default values defined in the function header.
