Emc <- function(dataset, popd.type = "sum", index = "saramaki", replicates=999, precision = 1e-9){
  #
  # The procedure will allow you to calculate likelihood measures of
  # niche breadth and overlap described in Petraitis (1979)
  #
  # Author: Nicola ZACCARELLI
  # E-mail: nicola.zaccarelli@gmail.com
  #
  # Version: 1.1
  # Date: 12/01/2018
  #
  # - Bar to show progress is removed.
  # - Introduced "precision": a filter to remove from the PS matrix all values which are
  #   lower than "precision". In this way some interactions are removed due to their small
  #   weight in the PS matrix.
  # - Introduced a way to handle isolated individuals in a network (thanks to Isabel's comment)
  #
 # for now only available for integer data
  # if (class(dataset) != "RInSp") stop("The input must be an object of class RInSp") # Changed because of the use of class(.)
  if (!inherits(dataset, "RInSp")) stop("The input must be an object of class RInSp")
  if (dataset$data.type != "integer") stop("Input data type must be integer.")
  if (popd.type %in% c("sum", "average") == FALSE) stop("The specified population diet type is wrong.")
  if (popd.type == "sum")
    { dietpop <- apply(dataset$resources, 2, sum) / sum(dataset$resources)} else {
      dietpop <- apply(dataset$proportion, 2, mean)}
  if (index %in% c("saramaki", "barrat") == FALSE) stop("The specified distance index is wrong.")
  if (index == "saramaki") d.index <- 1 else d.index <- 2
  if(!is.double(replicates)) replicates <- floor(abs(as.double(replicates)))
  if (replicates == 0) stop("The specified replicates value is wrong")
  if (replicates < 10)
      { replicates <- 10
        cat("\n Warning! Minimum number of replicates is 10. \n") }
  cat("\n Warning! resampling can take a while. Please be patient. \n")
totdieti <- apply(dataset$resources, 1, sum)
# pb <- txtProgressBar(min=0, max=10,  char= "+", style = 3)
#
# Change after Isabel comment
# We calculate the PS matrix and check for islated individuals
#
num.ind <- dataset$num.individuals
num.res <- dataset$num.prey
PS <- matrix(0, num.ind, num.ind)
for (i in 1:(num.ind - 1)) {
  for (j in (i+1):num.ind)
  {
    PS[i, j] <- 1 - 0.5 * sum(abs(dataset$proportions[i, ] - dataset$proportions[j, ]))
    if (PS[i,j] <= precision) PS[i,j] <- 0  # filter to remove small interactions
    PS[j, i] <- PS[i, j]
  }
}
if (sum(colSums(PS) == 0) > 0) Isolation <- 1 else Isolation <- 0
#
for (count in 1:10) {
     if (count < 10) StepCalc <- floor(replicates / 10) else StepCalc <- replicates - floor(replicates / 10) *9
     tmp <- .Call("CEmc", dataset$proportions, as.vector(d.index), as.vector(dietpop), as.vector(totdieti),
                 as.vector(StepCalc), as.vector(precision),PACKAGE="RInSp")
     if (count ==1) Ris <- tmp else Ris <- rbind(Ris, tmp[-1, ])
     # setTxtProgressBar(pb, count)
}
# close(pb)
attributes(Ris)$dimnames[[2]] <- c("Wmean", "E", "CW", "CwS")
cum.distr <- ecdf(Ris[, 2])
pvalue <- 1 - cum.distr(Ris[1, 2])
meannullE <- mean(Ris[-1, 2])
Eadj <- (Ris[1, 2] - meannullE) /(1 - meannullE)
CwS  <- Ris[1, 4]
cat("\n Araujo's E index (Eobs):", Ris[1,2])
cat("\n The mean Null E value is:", meannullE)
cat("\n The E adjusted value is:", Eadj)
cat("\n The p-value for P(Esim => Eobs) is:", pvalue)
cat("\n Degree of clustering in the network (Cws):", CwS)
# let's calculate the Monte Carlo resampling probability
cum.distr <- ecdf(Ris[, 4])
pvalue <- 1 - cum.distr(Ris[1, 4])
cat("\n The p-value for P(Cws.sim => Cws.obs) is:", pvalue)
cat("\n Population diet type :", popd.type)
cat("\n Cluster index :", index)
cat("\n")
if (Isolation > 0) cat("\n Please check your network because there are isolated elements.")
cat("\n")
Ris2 <- list(E = Ris[1, 2], meannullE= meannullE, Eadj= Eadj , p.value= pvalue, montecarlo= Ris, parameter = 2, pop.diet = popd.type,
             type.index = index, Precision = precision)
class(Ris2) <- "RInSp"
 return(Ris2)
}
