## crossing.R (2018-07-07)
##
## Functions for predicting genetic crosses between individuals
##
## Copyright 2018 Jan Engelstaedter
##
## This file is part of the R-package `peas'.


# Generates a table of all possible egg configurations
#
# In this table, each egg configuration is a way in which alleles could be drawn from a
# genotype to produce an egg haplotype. The table also has a column with expected
# fractions of these configurations.
# Given an actual genotype, each configuration can then be used to generate the actual egg haplotype.
#
# @param genopheno A genopheno object containing all information about the genetic setup.
#
# @return A data frame whose first columns represent all loci.
# For each locus, a number 1 represents the maternal allele, 2 represents the paternal allele,
# and 0 represents no allele. The last column contains the expected fraction of each configuration.
#
getEggConfigurations <- function(genopheno) {
  # generating table of possible configurations:
  poss <- list()  # list of possibilities for each locus
  for(lg in 1:length(genopheno$geno)) {
    if (genopheno$geno[[lg]]$type %in% c('autosomal', 'XY')) {
      poss <- c(poss, rep(list(c(1,2)), genopheno$geno[[lg]]$nloci))
    } else if (genopheno$geno[[lg]]$type == 'ZW') {
      poss <- c(poss, rep(list(c(1,0)), genopheno$geno[[lg]]$nloci))
    } else if (genopheno$geno[[lg]]$type == 'maternal') {
      poss <- c(poss, rep(list(c(1)), genopheno$geno[[lg]]$nloci))
    }
  }
  conf <- expand.grid(poss)

  # adding column names:
  confcolnames <- c()
  for(lg in 1:length(genopheno$geno)) {
    confcolnames <- c(confcolnames,
                      sapply(1:genopheno$geno[[lg]]$nloci, function(x) paste0(names(genopheno$geno)[lg],"_Loc",x) ))
  }
  colnames(conf) <- confcolnames

  # adding expected fraction for each configuration:
  conf$fraction <- NA
  for(i in 1:nrow(conf)) {
    frac <- 1 # fraction to be calculated
    ri <- 0  # current row index
    for(lg in 1:length(genopheno$geno)) {
      if (genopheno$geno[[lg]]$type %in% c('autosomal', 'XY', 'ZW')) {
        for(j in 1:genopheno$geno[[lg]]$nloci) {
          ri <- ri + 1
          if (j == 1) {
            frac <- frac * 0.5
          } else {
            if (conf[i, ri] != conf[i, ri - 1])
              frac <- frac * genopheno$geno[[lg]]$rec[j - 1]
            else frac <- frac * (1 - genopheno$geno[[lg]]$rec[j - 1])
          }
        }
      }
    }
    conf$fraction[i] <- frac
  }
  return(conf)
}

# Generates a table of all possible sperm configurations
#
# In this table, each sperm configuration is a way in which alleles could be drawn from a
# genotype to produce an egg haplotype. The table also has a column with expected
# fractions of these configurations.
# Given an actual genotype, each configuration can then be used to generate the actual sperm haplotype.
#
# @param genopheno A genopheno object containing all information about the genetic setup.
#
# @return A data frame whose first columns represent all loci.
# For each locus, a number 1 represents the maternal allele, 2 represents the paternal allele,
# and 0 represents no allele. The last column contains the expected fraction of each configuration.
#
getSpermConfigurations <- function(genopheno) {
  # generating table of possible configurations:
  poss <- list()  # list of possibilities for each locus
  for(lg in 1:length(genopheno$geno)) {
    if (genopheno$geno[[lg]]$type %in% c('autosomal', 'ZW')) {
      poss <- c(poss, rep(list(c(1,2)), genopheno$geno[[lg]]$nloci))
    } else if (genopheno$geno[[lg]]$type == 'XY') {
      poss <- c(poss, rep(list(c(1,0)), genopheno$geno[[lg]]$nloci))
    } else if (genopheno$geno[[lg]]$type == 'maternal') {
      poss <- c(poss, rep(list(c(0)), genopheno$geno[[lg]]$nloci))
    }
  }
  conf <- expand.grid(poss)

  # adding column names:
  confcolnames <- c()
  for(lg in 1:length(genopheno$geno)) {
    confcolnames <- c(confcolnames,
                      sapply(1:genopheno$geno[[lg]]$nloci, function(x) paste0(names(genopheno$geno)[lg],"_Loc",x) ))
  }
  colnames(conf) <- confcolnames

  # adding expected fraction for each configuration:
  conf$fraction <- NA
  for(i in 1:nrow(conf)) {
    frac <- 1 # fraction to be calculated
    ri <- 0  # current row index
    for(lg in 1:length(genopheno$geno)) {
      if (genopheno$geno[[lg]]$type %in% c('autosomal', 'XY', 'ZW')) {
        for(j in 1:genopheno$geno[[lg]]$nloci) {
          ri <- ri + 1
          if (j == 1) {
            frac <- frac * 0.5
          } else {
            if (conf[i, ri] != conf[i, ri - 1])
              frac <- frac * genopheno$geno[[lg]]$rec[j - 1]
            else frac <- frac * (1 - genopheno$geno[[lg]]$rec[j - 1])
          }
        }
      }
    }
    conf$fraction[i] <- frac
  }
  return(conf)
}


#' Predicts offspring distributions resulting from a cross between two genotypes.
#'
#' @param genopheno genopheno object specifying the genetic setup.
#' @param mom Genotype of the mother.
#' @param dad Genotype of the father.
#' @param equivalent Argument determining which genotypes are to be treated as equivalent.
#' Possible options are "phase" (genotypes that differ only by phase are equivalent,
#' i.e. the order of alleles within a locus does not matter),
#' "origin" (phase is important but not the origin of alleles),
#' and "none" (no two genotypes are equivalent).
#' @param output Specifies what the function should return: 'genotypes' returns the distribution of
#' offspring genotypes, 'phenotypes' returns the distribution of offspring phenotypes, and
#' 'both' returns both distributions.
#'
#' @return Depending on the output argument, a dataframe containing the distribution of
#' offspring genotypes, a dataframe containing the distribution of offspring phenotypes,
#' or both (combined in a list).
#' @export
#' @import stats
#'
#' @examples
#' MendelsPeas <- newGenopheno(alleleNames = list(c('Y', 'y')))
#' MendelsPeas <- setPhenotypes(MendelsPeas, 'colour', 'Y_', 'yellow')
#' MendelsPeas <- setPhenotypes(MendelsPeas, 'colour', 'yy', 'green')
#' predictCross(MendelsPeas, 'YY', 'yy')
#' predictCross(MendelsPeas, 'Yy', 'Yy')
#'
predictCross <- function(genopheno, mom, dad, equivalent = "phase", output = "both") {
  if (!isValidGenotype(mom, genopheno))
    stop("Maternal genotype not recognised.")
  if (!isValidGenotype(dad, genopheno))
    stop("Paternal genotype not recognised.")
  if (!(output %in% c("both","genotypes", "phenotypes")))
    stop("Unknown output argument.")
  if ((is.null(genopheno$pheno)) & (output != "genotypes")) {
    output <- "genotypes"
    warning("No phenotypes defined; only the offspring genotype distribution will be returned.")
  }

  momL <- convertGenoStringToList(mom, genopheno)  # maternal genotype in list format
  dadL <- convertGenoStringToList(dad, genopheno)  # paternal genotype in list format
  eggs <- getEggConfigurations(genopheno)
  sperm <- getSpermConfigurations(genopheno)
  noffspring <- nrow(eggs)*nrow(sperm)   # "raw" number of offspring genotypes
  if (is.null(genopheno$pheno)) { # no phenotypes defined
    offspring <- data.frame(matrix(NA, ncol=2, nrow=noffspring))
    colnames(offspring) <- c("fraction", "eqClass")
  } else {  # at least one trait defined
    offspring <- data.frame(matrix(NA, ncol=2 + ncol(genopheno$pheno), nrow=noffspring))
    colnames(offspring) <- c("fraction", "eqClass", colnames(genopheno$pheno))
  }
  offspringGTs <- vector("list", length = noffspring)
  io <- 0  # row in offspring dataframe

  # filling list of offspring genotypes and associated dataframe with fractions of these genotypes:
  for(ie in 1:nrow(eggs)) {
    for(is in 1:nrow(sperm)) {
      io <- io + 1  # index for offspring genotype
      offspring$fraction[io] <- eggs$fraction[ie] * sperm$fraction[is]
      offspringGTs[[io]] <- createEmptyGenotype(genopheno)
      j<-0
      # generating offspring genotype from egg and sperm configuration:
      for(lg in 1:length(genopheno$geno)) {  # linkage group in offspring genotype
        for(i in 1:genopheno$geno[[lg]]$nloci) {  # locus within linkage group in offspring genotype
          j <- j + 1
          offspringGTs[[io]][[lg]][[i]] <- c(momL[[lg]][[i]][eggs[ie,j]], dadL[[lg]][[i]][sperm[is,j]])
        }
      }
    }
  }

  # merging duplicates in the list and dataframe of offspring genotypes:
  ieq <- 0 # index for equivalence classes
  for(io in 1:noffspring) {
    if (is.na(offspring$eqClass[io])) {
      ieq <- ieq + 1
      equis <- getEquivalents(offspringGTs[[io]], genopheno, equivalent)
      offspring$eqClass[which(offspringGTs %in% equis)] <- ieq
    }
  }
  neq <- ieq # number of equivalency classes
  for(ieq in 1:neq) {
    firstMember <- match(ieq, offspring$eqClass)  # index of first genotype within equivalency class
    row.names(offspring)[firstMember] <- convertGenoListToString(offspringGTs[[firstMember]], genopheno)
    offspring[firstMember, "fraction"] <- sum(offspring[ieq == offspring$eqClass, "fraction"])
    if(!is.null(genopheno$pheno))
      offspring[firstMember, 3:ncol(offspring)] <- genopheno$pheno[row.names(offspring)[firstMember],]
  }
  offspring <- offspring[match(1:neq,offspring$eqClass), -2, drop = FALSE]

  # creating dataframe for phenotype distributions among offspring:
  if (output != "genotypes") {
    # replace all NA's first since they will otherwise be ignored:
    offspringPhenotypes <- offspring[-1]
    offspringPhenotypes[is.na(offspringPhenotypes)] <- "UnBeKnOwNsT"
    # group by distinct phenotype patterns:
    offspringPhenotypes <- aggregate(offspring[, "fraction", drop = FALSE],
                                     by = offspringPhenotypes,
                                     FUN = sum)
    # put NA's back in:
    offspringPhenotypes[offspringPhenotypes == "UnBeKnOwNsT"] <- NA
  }
  if (output == "both") {
    return(list(genotypes = offspring[1], phenotypes = offspringPhenotypes))
  } else if (output == "genotypes") {
    return(offspring[1])
  } else if (output == "phenotypes") {
    return(offspringPhenotypes)
  }
}
