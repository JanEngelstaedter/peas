## output.R (2018-07-07)
##
## Functions for displaying information about genetic systems
##
## Copyright 2018 Jan Engelstaedter
##
## This file is part of the R-package `peas'.


#' Print information about genotypes and corresponding phenotypes.
#'
#' @param x an obbject of class \code{genopheno}
#' @param ... further arguments passed to or from other methods
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
print.genopheno <- function(x, ...) {
  nLGs <- length(x$geno) # number of linkage groups
  if (nLGs == 1) {
    cat("Genetic system comprising one linkage group:\n")
  } else {
    cat("Genetic system comprising ",nLGs, " linkage groups:\n", sep="")
  }

  for(i in 1:nLGs) {
    if (x$geno[[i]]$nloci == 1) {
      cat("  Linkage group ",i,": ",x$geno[[i]]$type, ", ", x$geno[[i]]$nloci, " locus\n", sep="")
    } else if (x$geno[[i]]$nloci == 2) {
      cat("  Linkage group ",i,": ",x$geno[[i]]$type, ", ", x$geno[[i]]$nloci, " loci with recombination rate ", x$geno[[i]]$rec,"\n", sep="")
    } else if (x$geno[[i]]$nloci > 2) {
      cat("  Linkage group ",i,": ",x$geno[[i]]$type, ", ", x$geno[[i]]$nloci, " loci with recombination rates (", sep="")
      cat(x$geno[[i]]$rec, sep=", ")
      cat(")\n")
    }

    for(j in 1:x$geno[[i]]$nloci) {
      cat("    Alleles at locus ",j,": ", sep="")
      if (x$geno[[i]]$type %in% c('autosomal', 'maternal', 'paternal')) {
        for(k in 1:(x$geno[[i]]$nalleles[j] -1 ))
          cat(x$geno[[i]]$alleleNames[[j]][k],", ",sep="")
        cat(x$geno[[i]]$alleleNames[[j]][x$geno[[i]]$nalleles[j]],"\n",sep="")
      } else if (x$geno[[i]]$type %in% c('X-linked', 'Z-linked', 'Y-linked', 'W-linked')) {
        for(k in 1:(x$geno[[i]]$nalleles[j] -1 ))
          cat(x$geno[[i]]$alleleNames[[j]][k + 1],", ",sep="")
        cat(x$geno[[i]]$alleleNames[[j]][x$geno[[i]]$nalleles[j] + 1],"\n",sep="")
      }
    }
  }

  if (is.null(x$pheno)) {
    cat('No phenotypes defined.')
  } else {
    cat('Phenotypes defined for the following traits:\n')
    for(i in 1:ncol(x$pheno)) {
      cat("  ", colnames(x$pheno)[i], " (trait values: ", sep="")
      cat(unique(x$pheno[,i]), sep=", ")
      cat(")\n")
    }
    if(any(is.na(x$pheno))) cat("\n  Note: some genotypes have undefined trait values.")
  }
}


#' Retrieve information about the phenotypes of the genetic system
#'
#' @param genopheno a genopheno object containing all information about a genetic system.
#' @param equivalent argument determining which genotypes are to be treated as equivalent.
#' Possible options are "phase" (genotypes that differ only by phase are equivalent,
#' i.e. the order of alleles within a locus does not matter),
#' "origin" (phase is important but not the origin of alleles),
#' and "none" (no two genotypes are equivalent).
#' @param hideNoncanonical if this argument is set to TRUE (the default value),
#' genotypes involving sex chromomal linkage groups that normally cannot arise
#' (e.g., complete lack of an allele at an X-linked locus) are not displayed.
#' If this argument is set to FALSE, even these noncanonical genotypes will be displayed.
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
getPhenotypes <- function(genopheno, equivalent = "phase", hideNoncanonical = TRUE) {
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
        # suppressing non-canonical genotypes in case of X- or Z-linked inheritance:
        if (hideNoncanonical & (isNoncanonical(convertGenoStringToList(rownames(genopheno$pheno)[i], genopheno), genopheno)))
          toShow[i] <- FALSE
      }
      for(i in 1:nrow(genopheno$pheno)) {
        if (toShow[i]) {
          # suppressing equivalent genotypes:
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
