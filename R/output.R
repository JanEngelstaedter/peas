

#' Print information about genotypes and corresponding phenotypes.
#'
#' @param gp Object of class \code{genopheno}
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' MendelsPeas <- newGenopheno(nloci = 1, alleleNames = list(c('Y', 'y')))
#' MendelsPeas <- setPhenotypes(MendelsPeas, 'colour', 'Y_', 'yellow')
#' MendelsPeas <- setPhenotypes(MendelsPeas, 'colour', 'yy', 'green')
#' MendelsPeas
#'
print.genopheno <- function(gp) {
  nLGs <- length(gp$geno) # number of linkage groups
  if (nLGs == 1)
    cat("Genetic system comprising one linkage group:\n")
  else
    cat("Genetic system comprising ",nLGs, " linkage groups:\n", sep="")
  for(i in 1:nLGs) {
    if (gp$geno[[i]]$nloci == 1)
      cat("  Linkage group ",i,": ",gp$geno[[i]]$type, ", ", gp$geno[[i]]$nloci, " locus\n", sep="")
    else if (gp$geno[[i]]$nloci == 2)
      cat("  Linkage group ",i,": ",gp$geno[[i]]$type, ", ", gp$geno[[i]]$nloci, " loci with recombination rate ", gp$geno[[i]]$rec,"\n", sep="")
    else if (gp$geno[[i]]$nloci > 2) {
      cat("  Linkage group ",i,": ",gp$geno[[i]]$type, ", ", gp$geno[[i]]$nloci, " loci with recombination rates (", sep="")
      cat(gp$geno[[i]]$rec, sep=", ")
      cat(")\n")
    }

    for(j in 1:gp$geno[[i]]$nloci) {
      cat("    Alleles at locus ",j,": ", sep="")
      for(k in 1:(gp$geno[[i]]$nalleles[j] -1 ))
        cat(gp$geno[[i]]$alleleNames[[j]][k],", ",sep="")
      cat(gp$geno[[i]]$alleleNames[[j]][gp$geno[[i]]$nalleles[j]],"\n",sep="")
    }
  }

  if (is.null(gp$pheno)) {
    cat('No phenotypes defined.')
  } else {
    cat('Phenotypes defined for the following traits:\n')
    for(i in 1:ncol(gp$pheno)) {
      cat("  ", colnames(gp$pheno)[i], " (trait values: ", sep="")
      cat(unique(gp$pheno[,i]), sep=", ")
      cat(")\n")
    }
    if(any(is.na(gp$pheno))) cat("\n  Note: some genotypes have undefined trait values.")
  }
}


#' Retrieve information about the phenotypes of the genetic system
#'
#' @param genopheno A genopheno object containing all information about a genetic system.
#' @param equivalent Argument determining which genotypes are to be treated as equivalent.
#' Possible options are "phase" (genotypes that differ only by phase are equivalent,
#' i.e. the order of alleles within a locus does not matter),
#' "origin" (phase is important but not the origin of alleles),
#' and "none" (no two genotypes are equivalent).
#'
#' @return A dataframe in which rows represent genotypes (given by row names),
#'         and columns represent traits.
#' @export
#'
#' @examples
#'
#' MendelsPeas <- newGenopheno(nloci = 1, alleleNames = list(c('Y', 'y')))
#' MendelsPeas <- setPhenotypes(MendelsPeas, 'colour', 'Y_', 'yellow')
#' MendelsPeas <- setPhenotypes(MendelsPeas, 'colour', 'yy', 'green')
#' getPhenotypes(MendelsPeas)
#'
getPhenotypes <- function(genopheno, equivalent = "phase") {
  if (is.null(genopheno$pheno)) {
    print("No phenotypes specified yet.")
  } else {
    if (equivalent == 'none') {
      return(genopheno$pheno)
    } else if ((equivalent == 'phase') | (equivalent == 'origin')) {
      issueWarning <- FALSE
      toShow <- rep(TRUE, nrow(genopheno$pheno))
      names(toShow) <- rownames(genopheno$pheno)
      for(i in 1:nrow(genopheno$pheno)) {
        if (toShow[i]) {
          equisList <- getEquivalents(convertGenoStringToList(rownames(genopheno$pheno)[i], genopheno),
                                      genopheno, equivalent = equivalent) # list of genotypes that are equivalent to gt #i
          equisString <- sapply(equisList, convertGenoListToString, genopheno = genopheno) # same in string format
          for(j in 1:ncol(genopheno$pheno))
            if (length(unique(genopheno$pheno[equisString,j])) > 1) issueWarning <- TRUE
          equisString <- equisString[equisString!=rownames(genopheno$pheno)[i]]   # remove current gt from vector
          toShow[equisString] <- FALSE
        }
      }
      if (issueWarning) warning("Genotypes considered equivalent have been assigned different phenotypes.")
      return(genopheno$pheno[toShow,,drop = FALSE])
    } else {
      stop("Unknown value for argument 'equivalent")
    }
  }
}
