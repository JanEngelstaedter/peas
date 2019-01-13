## setup.R (2018-07-07)
##
## Functions for defining genetic systems
##
## Copyright 2018 Jan Engelstaedter
##
## This file is part of the R-package `peas'.

#' Generates a new genetic system
#'
#' This function generates a new genetic system that stores information about the number of loci, alleles at each locus
#' as well as recombination rates and mode of inheritance.
#'
#' @param nloci The number of loci. The default value is 1 but any (reasonable) number can be chosen.
#' @param nalleles The number of alleles for each locus. This should be a vector of length \code{nloci} containing positive integer numbers.
#' By default, all loci are set to be biallelic.
#' @param alleleNames The names for the alleles at the different loci.
#' The default value is 'standard', which for biallelic loci produces allele names 'a' and 'A' for the first locus, 'b' and 'B' for the second locus etc.
#' If one locus has more than two alleles, 'standard' produces allele names a0, a1, a2; b0, b1, b2 etc.
#' For other allele names, \code{alleleNames} should be a list of nloci elements, each of which needs to be a vector containing the allele names.
#' Allele names need to be strings (class \code{character}) and can contain letters, numbers,
#' or symbols (except '~', '.' and '|').
#' @param rec A vector of length \code{nloci-1} specifying recombination rates between the loci.
#' @param type This parameter specifies the ploidy and inheritance for the loci.
#' The default value is "autosomal" and currently, this is also the only supported type
#' of linkage group.
#'
#' @return An object of class \code{genopheno} that stores all of the specifications provided about the genetic system
#' and that can then be extended to also include other linkage groups and corresponding phenotypes.
#' @export
#'
#' @examples
#' MendelsPeas <- newGenopheno(nloci = 1, alleleNames = list(c('Y', 'y')))
#' MendelsPeas <- addLinkageGroup(MendelsPeas, nloci = 1, alleleNames = list(c('R', 'r')))
#' MendelsPeas <- setPhenotypes(MendelsPeas, 'colour', 'Y_ | __', 'yellow')
#' MendelsPeas <- setPhenotypes(MendelsPeas, 'colour', 'yy | __', 'green')
#' MendelsPeas <- setPhenotypes(MendelsPeas, 'shape', '__ | R_', 'round')
#' MendelsPeas <- setPhenotypes(MendelsPeas, 'shape', '__ | rr', 'wrinkled')
#' MendelsPeas
#'
newGenopheno <- function(nloci = 1, nalleles = rep(2, nloci), alleleNames = 'standard', rec = NULL, type = 'autosomal') {
  # initial checking if all arguments make sense and are consistent:
  checkGenophenoArguments(nloci, nalleles, alleleNames, rec, type)
  gp <- list(geno = list(), pheno = NULL)
  class(gp) <- "genopheno"
  gp <- addLinkageGroup(gp, nloci, nalleles, alleleNames, rec, type)
  return(gp)
}


#' Add linkage group to genetic system
#'
#' @param genopheno genopheno object specifying the genetic setup.
#' @inheritParams newGenopheno
#'
#' @return An object of class \code{genopheno} that stores all of the specifications provided about the genetic setup
#' and that can then be extended to also include other linkage groups and corresponding phenotypes.
#' @export
#'
addLinkageGroup <- function(genopheno, nloci = 1, nalleles = rep(2, nloci), alleleNames = 'standard', rec = NULL, type = 'autosomal') {
  checkGenophenoArguments(nloci, nalleles, alleleNames, rec, type)
  if (!is.null(genopheno[["pheno"]]))
    stop("New linkage groups can only be added before phenotypes are specified.", call. = FALSE)
  nLGs <- length(genopheno[["geno"]])
  if (alleleNames[1]=='standard') alleleNames <- getStandardAlleleNames(nloci, nalleles)
  if (type %in% c('X-linked', 'Z-linked', 'Y-linked', 'W-linked')) {
    # check if there is already a sex chromosome linkage group:
    if ((nLGs > 0) & (length(intersect(sapply(genopheno$geno, with, type), c('X-linked', 'Z-linked', 'Y-linked', 'W-linked'))) > 0))
      stop("Only one sex chromosome linkage group can be defined.", call. = FALSE)

    for(i in 1:nloci)
      alleleNames[[i]] <- c(strrep(".", nchar(alleleNames[[i]][1])), alleleNames[[i]])
  }
  genopheno[["geno"]][[nLGs + 1]] <- list(nloci = as.integer(nloci),
                                          nalleles = as.integer(nalleles),
                                          alleleNames = alleleNames,
                                          rec = rec,
                                          type = type)
  names(genopheno[["geno"]])[nLGs + 1] <- paste0("LG",nLGs + 1)
  return(genopheno)
}

# This function performs some checks to see if all arguments supplied to the functions
# newGenopheno and addLinkageGroup are valid.
checkGenophenoArguments <- function(nloci, nalleles, alleleNames, rec, type) {
  if ((!is.numeric(nloci)) | (length(nloci)>1) | (nloci < 1))
    stop("Number of loci must be a positive integer number.", call. = FALSE)
  if (length(nalleles) != nloci)
    stop(paste0("Vector of length nloci=",nloci," expected for argument 'nalleles'."), call. = FALSE)
  if (alleleNames[1] != 'standard') {
    if ((!is.list(alleleNames)) | (length(alleleNames)!=nloci))
      stop(paste0("List of length nloci=", nloci, " expected for argument 'alleleNames'."), call. = FALSE)
    if (!identical(sapply(alleleNames, length), as.integer(nalleles)))
      stop("Number of alleles does not correspond to allele names.", call. = FALSE)
    if (any(sapply(alleleNames, function(x) { length(unique(nchar(x)))>1 })))
      stop("Allele names need to have the same number of characters for each locus.", call. = FALSE)
    if (any(sapply(alleleNames, function(x) { any(unique(x) != x) })))
      stop("For each locus allele names need to be all different.", call. = FALSE)
    if (any(sapply(alleleNames, function(x) { any(sapply(c('_','~','\\|','\\.',' '), grepl, x = x)) })))
      stop("Invalid character in one of the allele names.", call. = FALSE)
  }
  if ((nloci == 1) & (!is.null(rec)))
    stop("No recombination rates can be set when there is only one locus.", call. = FALSE)
  if ((type %in% c('maternal', 'paternal')) & (!is.null(rec)))
    stop("No recombination rates can be set for uniparentally inherited linkage groups.", call. = FALSE)
  if ((type %in% c('autosomal', 'X-linked', 'Z-linked', 'Y-linked', 'W-linked')) & (!is.null(rec)) & (length(rec) != (nloci-1)))
    stop(paste("Vector of lenght",nloci-1,"expected for argument 'rec'."), call. = FALSE)
  if (!(type %in% c('autosomal', 'X-linked', 'Z-linked', 'Y-linked', 'W-linked', 'maternal', 'paternal')))
    stop("Inheritance specified by 'type' not supported.", call. = FALSE)
}



# Produces a list of standard allele names for a locus
#
# @inheritParams newGenopheno
#
# @return A vector of character strings.
#
getStandardAlleleNames <- function(nloci, nalleles) {
  alleleNames<-vector("list", nloci)
  if (max(nalleles)==2) {
    for (locus in 1:nloci) alleleNames[[locus]]<-c(letters[locus],LETTERS[locus])
  } else if (max(nalleles)>2 & max(nalleles)<=10) {
    for (locus in 1:nloci) alleleNames[[locus]]<-paste0(letters[locus],0:(nalleles[locus] -1))
  }
  else stop("Standard allele names not available when there are more than ten alleles at a locus.")
  return(alleleNames)
}


#' Set trait values for one or several genotypes.
#'
#' @param genopheno genopheno object specifying the genetic setup.
#' @param genotypes A vector of strings specifying genotypes.
#' @param traitName The name of the trait whose values are to be set.
#' @param traitValue The trait value.
#' @param equivalent Argument determining which genotypes are to be treated as equivalent.
#' Possible options are "phase" (genotypes that differ only by phase are equivalent,
#' i.e. the order of alleles within a locus does not matter),
#' "origin" (phase is important but not the origin of alleles),
#' and "none" (no two genotypes are equivalent).
#'
#' @return genopheno object with trait information added.
#' @export
#'
#' @examples
#' MendelsPeas <- newGenopheno(nloci = 1, alleleNames = list(c('Y', 'y')))
#' MendelsPeas <- setPhenotypes(MendelsPeas, 'colour', 'Y_', 'yellow')
#' MendelsPeas <- setPhenotypes(MendelsPeas, 'colour', 'yy', 'green')
#' getPhenotypes(MendelsPeas)
#'
setPhenotypes <- function(genopheno, traitName, genotypes, traitValue, equivalent = "phase") {
  genotypes <- tidyGenoString(genotypes)
  if (is.null(genopheno$pheno)) {   # first time phenotypes are specified
    allGenoStrings <- getAllGenoStrings(genopheno)
    genopheno$pheno <- data.frame(rep(NA, length(allGenoStrings)), row.names = allGenoStrings)
    colnames(genopheno$pheno)[1] <- traitName
  } else {  # some phenotypes have already been specified
    if (!(traitName %in% colnames(genopheno$pheno)))  # new trait
      genopheno$pheno[,traitName] <- NA
  }

  for(i in 1:length(genotypes)) {
    phs<-unlist(gregexpr(pattern='_',genotypes[i])) # placeholder positions
    if (phs[1]>=1)  # are there any placeholders?
    {
      gwp<-insertPlaceholders(rownames(genopheno$pheno), phs)  # all genotypes with placeholders inserted
      matches<-rownames(genopheno$pheno)[which(genotypes[i]==gwp)]
    } else {
      if (genotypes[i] %in% rownames(genopheno$pheno)) matches<-genotypes[i]
    }

    if(length(matches)>0)
    {
      for(j in 1:length(matches)) {
        equis <- sapply(getEquivalents(convertGenoStringToList(matches[j], genopheno), genopheno, equivalent = equivalent),
                        convertGenoListToString, genopheno = genopheno)
        genopheno$pheno[equis, traitName]<-traitValue
      }
    } else warning(paste0("No matching genotypes found for ", genotypes[i], "."))
  }
  return(genopheno)
}


#' Add sex phenotypes to a genetic system with sex chromosomes
#'
#' This function add a trait 'sex' to an existing genopheno object.
#' Trait values ('female' and 'male') are assigned automatically to all genotypes
#' according to sex chromosome configurations.
#'
#' @param genopheno genopheno object specifying the genetic setup.
#'
#' @return genopheno object with sex trait information added.
#' @export
#'
#' @examples
#' mySystem <- newGenopheno(type = "X-linked")
#' mySystem <- setSexPhenotypes(mySystem)
#' getPhenotypes(mySystem)
#'
setSexPhenotypes <- function(genopheno) {
  # identify sex chromosomes:
  sexCSystem <- intersect(sapply(genopheno$geno, with, type), c('X-linked', 'Z-linked', 'Y-linked', 'W-linked'))
  if (length(sexCSystem) == 0) # no sex chromosomes
    stop("No sex chromosome linkage groups defined.")
  if (length(sexCSystem) > 1) # more than one sex chromosome system
    stop("Sex chromosome system ambiguous.")
  maleGTs <- ""
  femaleGTs <- ""
  for(lg in 1:length(genopheno$geno)) {
    if (lg > 1) {
      maleGTs <- paste0(maleGTs, "|")
      femaleGTs <- paste0(femaleGTs, "|")
    }
    if (genopheno$geno[[lg]]$type == 'autosomal')
      for(i in 1:genopheno$geno[[lg]]$nloci) {
        if (i>1) {
          maleGTs <- paste0(maleGTs, "~")
          femaleGTs <- paste0(femaleGTs, "~")
        }
        maleGTs <- paste0(maleGTs, paste0(rep("_", 2*nchar(genopheno$geno[[lg]]$alleleNames[[i]][1])), collapse = ""))
        femaleGTs <- paste0(femaleGTs, paste0(rep("_", 2*nchar(genopheno$geno[[lg]]$alleleNames[[i]][1])), collapse = ""))
      }
    if (genopheno$geno[[lg]]$type == 'X-linked')
      for(i in 1:genopheno$geno[[lg]]$nloci) {
        if (i>1) {
          maleGTs <- paste0(maleGTs, "~")
          femaleGTs <- paste0(femaleGTs, "~")
        }
        maleGTs <- paste0(maleGTs, paste0(rep("_", nchar(genopheno$geno[[lg]]$alleleNames[[i]][1])), collapse = ""))
        maleGTs <- paste0(maleGTs, paste0(rep(".", nchar(genopheno$geno[[lg]]$alleleNames[[i]][1])), collapse = ""))
        femaleGTs <- paste0(femaleGTs, paste0(rep("_", 2*nchar(genopheno$geno[[lg]]$alleleNames[[i]][1])), collapse = ""))
      }
    if (genopheno$geno[[lg]]$type == 'Z-linked')
      for(i in 1:genopheno$geno[[lg]]$nloci) {
        if (i>1) {
          maleGTs <- paste0(maleGTs, "~")
          femaleGTs <- paste0(femaleGTs, "~")
        }
        maleGTs <- paste0(maleGTs, paste0(rep("_", 2*nchar(genopheno$geno[[lg]]$alleleNames[[i]][1])), collapse = ""))
        femaleGTs <- paste0(femaleGTs, paste0(rep(".", nchar(genopheno$geno[[lg]]$alleleNames[[i]][1])), collapse = ""))
        femaleGTs <- paste0(femaleGTs, paste0(rep("_", nchar(genopheno$geno[[lg]]$alleleNames[[i]][1])), collapse = ""))
      }
    if (genopheno$geno[[lg]]$type == 'Y-linked')
      for(i in 1:genopheno$geno[[lg]]$nloci) {
        if (i>1) {
          maleGTs <- paste0(maleGTs, "~")
          femaleGTs <- paste0(femaleGTs, "~")
        }
        maleGTs <- paste0(maleGTs, paste0(rep(".", nchar(genopheno$geno[[lg]]$alleleNames[[i]][1])), collapse = ""))
        maleGTs <- paste0(maleGTs, paste0(rep("_", nchar(genopheno$geno[[lg]]$alleleNames[[i]][1])), collapse = ""))
        femaleGTs <- paste0(femaleGTs, paste0(rep(".", 2*nchar(genopheno$geno[[lg]]$alleleNames[[i]][1])), collapse = ""))
      }
    if (genopheno$geno[[lg]]$type == 'W-linked')
      for(i in 1:genopheno$geno[[lg]]$nloci) {
        if (i>1) {
          maleGTs <- paste0(maleGTs, "~")
          femaleGTs <- paste0(femaleGTs, "~")
        }
        maleGTs <- paste0(maleGTs, paste0(rep(".", 2*nchar(genopheno$geno[[lg]]$alleleNames[[i]][1])), collapse = ""))
        femaleGTs <- paste0(femaleGTs, paste0(rep("_", nchar(genopheno$geno[[lg]]$alleleNames[[i]][1])), collapse = ""))
        femaleGTs <- paste0(femaleGTs, paste0(rep(".", nchar(genopheno$geno[[lg]]$alleleNames[[i]][1])), collapse = ""))
      }
  }
  if (sexCSystem %in% c('X-linked', 'W-linked')) {
    genopheno <- setPhenotypes(genopheno, "sex", femaleGTs, "female")
    genopheno <- setPhenotypes(genopheno, "sex", maleGTs, "male")
  } else {
    genopheno <- setPhenotypes(genopheno, "sex", maleGTs, "male")
    genopheno <- setPhenotypes(genopheno, "sex", femaleGTs, "female")
  }
  return(genopheno)
}
