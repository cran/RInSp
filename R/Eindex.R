Eindex <-
  function(dataset, index = "saramaki", jackknife = FALSE, precision = 1e-9){
    #
    # The procedure calculates the E measure of interindividual variation,
    # its variance and the value of the Cws measure of modularity following Araujo et al. (2008).
    #
    # Author: Nicola ZACCARELLI
    # E-mail: nicola.zaccarelli@gmail.com
    #
    # Version: 1.2
    # Date: 12/01/2018
    #
    # - Bar to show progress is removed.
    # - Introduced "precision": a filter to remove from the PS matrix all values which are
    #   lower than "precision". In this way some interactions are removed due to their small
    #   weight in the PS matrix.
    # - Introduced a way to handle isolated individuals in a network (thanks to Isabel's comment)
    #
    # some checking
    # if (class(dataset) != "RInSp") stop("The input must be an object of class RInSp") # Changed because of the use of class(.)
    if (!inherits(dataset, "RInSp")) stop("The input must be an object of class RInSp")
    if ((index != "saramaki") & (index != "barrat")) stop("Wrong index type. Must be Saramaki or Barrat.")
    if (jackknife %in% c("T", "F", "TRUE", "FALSE") == FALSE) stop("Wrong jackknife option. Must be TRUE or FALSE.")

    # Calculate values for O, mean O and E
    num.ind <- dataset$num.individuals
    num.res <- dataset$num.prey
    PS <- matrix(0, num.ind, num.ind)
    O <- 0
    for (i in 1:(num.ind - 1)) {
      for (j in (i+1):num.ind)
      {
        PS[i, j] <- 1 - 0.5 * sum(abs(dataset$proportions[i, ] - dataset$proportions[j, ]))
        if (PS[i,j] <= precision) PS[i,j] <- 0  # filter to remove small interactions
        PS[j, i] <- PS[i, j]
        O <- PS[i, j] + O
      }
    }
    Oaverage <- 2 * O / ((num.ind - 1)* num.ind)

    # Calculate variance based on jackknife (i.e., a leave one individual out approach)
    if (jackknife == TRUE) {
      # cat("\n Jackknife in progress... \n")
      # # Progress bar with 10% increment intervals
      # pb <- txtProgressBar(min=0, max=10,  char= "+", style <- 3)
      # stepwidth <- floor(dataset$num.individuals / 10)
      # step <- 1

      num.ind2 <- num.ind - 1
      num.res2 <- num.res - 1

      for (i in 1:dataset$num.individuals)
      {
        # if (i== step * stepwidth) {
        #   setTxtProgressBar(pb, step)
        #   step <- step + 1
        # }
        indkeep <- c(1:dataset$num.individuals) != i
        datatmppro <- subset(dataset$proportions, subset=indkeep)
        PS2 <- matrix(0, num.ind2, num.ind2)
        O2 <- 0
        for (k in 1:(num.ind2 - 1)) {
          for (j in (k+1):num.ind2){
            PS2[k, j] <- 1 - 0.5 * sum(abs(datatmppro[k, ] - datatmppro[j, ]))
            if (PS2[k,j] <= precision) PS2[k,j] <- 0  # filter to remove small interactions
            PS2[j, k] <- PS2[k, j]
            O2 <- PS2[k, j] + O2
          }
        }
        if (i == 1) Ejack <- O2 else Ejack <- c(Ejack, O2)
      }

      Ejack <- 1 - 2 * Ejack/ ((num.ind2 - 1)* num.ind2)
      vare <- (sum(Ejack^2) - ((sum(Ejack))^2 / num.ind)) / (num.ind * (num.ind - 1))
      # close(pb)
    } else {
      Ejack <- 0
      vare <- 0
    }

    # cat("\n Please wait, we are computing E and Cws. \n")

    #Calculate the cluster index
    Wmax <- max(PS)
    # # Progress bar with 10% increment intervals
    # pb <- txtProgressBar(min=0, max=10,  char= "+", style <- 3)
    # stepwidth <- floor(dataset$num.individuals / 10)
    # step <- 1

    if (index == "saramaki")
    {
      PSbinary <- matrix(0, num.ind, num.ind)
      PSbinary[PS > 0] <- 1
      Ki <- apply(PSbinary, 1, sum)
      for (i in 1:num.ind)
       { # if (i== step * stepwidth) {
      #   setTxtProgressBar(pb, step)
      #   step <- step + 1
      # }
      F=0
      for (j in 1:num.ind)
      {
        for (k in 1:num.ind)
        {
          # When in a clique there is a zero link element, the product is zero. For this those two line are not necessary
          #           if ((i != j) && (i != k) && (j != k))
          #           {
          #             if ((PS[i, j] > 0) && (PS[i, k] > 0) && (PS[j, k] > 0))
          #             {
          Wraiz3 <- ((PS[i, j]*PS[j, k]*PS[k, i])^(1/3))/Wmax
          F <- F + Wraiz3;
          #              }
          #             }
        }
      }
      if  (Ki[i] <= 1)
      { temp <- 1                             # Georgia change
        Ki[i] <- 1 }
      else temp <- (Ki[i] - 1.0) * Ki[i]
      if (i == 1) Cw <- F / temp else Cw <- c(Cw, F / temp)
      # Change prompted by Isabel to deal with isolated nodes
      if (sum(Cw == 0) > 0) Isolation <- 1 else Isolation <- 0
      } # For
    } else  # barrat index
    {
      PSbinary <- matrix(0, num.ind, num.ind)
      PSbinary[PS > 0] <- 1
      Ki <- apply(PSbinary, 1, sum)
      PSbinary <- PS
      PSbinary[PS <= 0] <- 0
      Si <- apply(PSbinary, 1, sum)
      for (i in 1:num.ind)
      { # if (i== step * stepwidth) {
        # setTxtProgressBar(pb, step)
        # step <- step + 1
      # }
      F=0;
      for (j in 1:num.ind)
      {
        for (k in 1:num.ind)
        {
           if ((i != j) && (i != k) && (j != k))
          {
            if ((PS[i, j] > 0) && (PS[i, k] > 0) && (PS[j, k] > 0))
              {
               F <- F + (PS[i, j] + PS[i, k])/2.0
              }
           }
        }
      }
      temp <- (Ki[i] - 1.0)
      if (temp <= 0) temp <- 1.0
      if (i == 1) Cw <- F /(Si[i]*temp) else Cw <- c(Cw, F /(Si[i]*temp))
      }
    }
   # close(pb)
    # Change prompted by Isabel to deal with isolated nodes
    if (sum(!is.finite(Cw)) > 0) {
      Isolation <- 1
      Cw[!is.finite(Cw)] <- 0
    }  else Isolation <- 0
    CW <- sum(Cw)/num.ind;
    if ((CW + Oaverage) != 0) CwS <- (CW - Oaverage)/(CW + Oaverage)

    # Calculate the binary matrix
    PSbinary <- matrix(0, num.ind, num.ind)
    PSbinary[PS >= O] <- 1
    Ris <- list(Omean = Oaverage, E = 1 - Oaverage, PS = PS, PSbinary = PSbinary, Ejack = Ejack, VarE = vare,
               CW = CW, CwS = CwS, Cw= Cw, index = index, Ki = Ki, Precision = precision, Isolation = Isolation)
    cat("\n Araujo's E:", 1 - Oaverage)
    cat("\n F value:", sprintf("%.10f", F))
#    cat("\n Ki", sprintf("%.1f", Ki))
#    cat("\n Cws ", sprintf("%.10f", Cw))
    if (jackknife == TRUE) cat("\n Araujo's E variance:", vare)
    cat("\n Degree of clustering in the network (Cws):", CwS)
    cat("\n Index type:", index)
    if (Isolation > 0) cat("\n Please check your network because there are isolated elements.")
    cat("\n")
    return(Ris)
  }
