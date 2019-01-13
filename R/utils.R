## utils.R (2018-07-07)
##
## Various non-exported helper functions
##
## Copyright 2018 Jan Engelstaedter
##
## This file is part of the R-package `peas'.


#' Genotype string tidying
#'
#' Converts a character string representing a genotype into a canonical form
#' where spaces within linkage groups are removed and linkage groups are
#' connected by ' | '.
#'
#' @param genoString The character string to be tidied.
#'
#' @return A tidy genotype string.
#'
#' @keywords internal
#'
#' tidyGenoString("  aA ~bb|CC")
tidyGenoString <- function(genoString) {
  tidyGeno <- gsub(" ", "", genoString)
  tidyGeno <- gsub("\\|", " | ", tidyGeno)
  tidyGeno
}


# insert placeholders ("_") in all strings at positions
# specified by the vector pos
insertPlaceholders<-function(strings, pos)
{
  for(p in pos) substr(strings,p,p) <- "_"
  return(strings)
}


#' Convert from genotype string to list format
#'
#' Genotypes are internally represented as lists but specified as strings by the user.
#' This function converts from the latter into the former format.
#'
#' @param genoString Character string representing a genotype
#' @param genopheno The genetic system that the genotype is part of
#'
#' @keywords internal
#'
#' @return A genotype object in list format
#'
convertGenoStringToList <- function(genoString, genopheno) {
  genoString <- tidyGenoString(genoString)  # removing excess spaces etc.
  nLGs <- length(genopheno$geno)  # expected number of linkage groups
  genoList <- vector("list", nLGs)
  spos <- 1  # position within genoString
  for(lg in 1:nLGs) {  # loop over all linkage groups
    if (genopheno$geno[[lg]]$type %in% c('autosomal', 'X-linked', 'Z-linked', 'Y-linked', 'W-linked')) {
      ploidy <- 2
    } else if (genopheno$geno[[lg]]$type %in% c('maternal', 'paternal')) {
      ploidy <- 1
    } else {
      stop("Linkage groups not supported")
    }
    if (lg > 1) spos <- spos + 3 # skipping the " | " seperator
    genoList[[lg]] <- vector("list", genopheno$geno[[lg]]$nloci)
    for(i in 1:genopheno$geno[[lg]]$nloci) { # loop over all loci
      len <- nchar(genopheno$geno[[lg]]$alleleNames[[i]][1])
      if (i>1) spos <- spos + 1  # skip tilde
      genoList[[lg]][[i]] <- rep(NA, ploidy)
      for(pos in 1:ploidy) { # loop over the alleles per locus
        for(j in 1:length(genopheno$geno[[lg]]$alleleNames[[i]])) { # loop over allele names
          if (substr(genoString, spos, spos + len -1) == genopheno$geno[[lg]]$alleleNames[[i]][j]) {  # match!
            genoList[[lg]][[i]][pos] <- j
            if (genopheno$geno[[lg]]$type %in% c('X-linked', 'Z-linked', 'Y-linked', 'W-linked'))
              genoList[[lg]][[i]][pos] <- genoList[[lg]][[i]][pos] - 1  # so that '.' becomes assigned a zero
            break
          }
        }
        if(is.na(genoList[[lg]][[i]][pos]))
          stop(paste0("Allele in linkage group ", lg, ", locus ",i," not recognised."))
        spos <- spos + len
      }
    }
  }
  return(genoList)
}

#' Convert from genotype list to string format
#'
#' Genotypes are internally represented as lists but specified as strings by the user.
#' This function converts from the former into the latter format.
#'
#' @param genoList A genotype representation in list format.
#' @param genopheno The genetic system that the genotype is part of
#'
#' @keywords internal
#'
#' @return Genotype in string format
#'
convertGenoListToString <- function(genoList, genopheno) {
  genoString <- "" # start with empty string
  for(lg in 1:length(genoList)) { # loop over all linkage groups
    if(lg>1) genoString <- paste0(genoString, " | ")
    if (genopheno$geno[[lg]]$type %in% c('autosomal', 'maternal', 'paternal')) {
      for(i in 1:length(genoList[[lg]])) { # loop over all loci
        if(i>1) genoString <- paste0(genoString, "~")
        for(j in 1:length(genoList[[lg]][[i]]))  # loop over alleles at locus
          genoString <- paste0(genoString, genopheno$geno[[lg]]$alleleNames[[i]][[ genoList[[lg]][[i]][j] ]])
      }
    } else if (genopheno$geno[[lg]]$type %in% c('X-linked', 'Z-linked', 'Y-linked', 'W-linked')) {
      for(i in 1:length(genoList[[lg]])) { # loop over all loci
        if(i>1) genoString <- paste0(genoString, "~")
        for(j in 1:length(genoList[[lg]][[i]]))  # loop over alleles at locus
          genoString <- paste0(genoString, genopheno$geno[[lg]]$alleleNames[[i]][[ genoList[[lg]][[i]][j] + 1 ]])
      }
    }
  }
  return(genoString)
}

# Swapping alleles at the two positions of a locus
#
# @param gt Genotype in list format
# @param lg Linkage group numer
# @param loc Locus number
#
# @return Genotype in list format
#
swapAlleles <- function(gt, lg, loc) {
  dummy <- gt[[lg]][[loc]][1]
  gt[[lg]][[loc]][1] <- gt[[lg]][[loc]][2]
  gt[[lg]][[loc]][2] <- dummy
  return(gt)
}


# Swapping the two alleles at all loci of a linkage group
#
# @param gt Genotype in list format.
# @param lg Linkage group where alleles should be swapped.
#
# @return A genotype in list format.
#         This genotype is identical to the supplied one except for the swapped linkage group.
#
swapAllAlleles <- function(gt, lg) {
  for(i in 1:length(gt[[lg]])) gt <- swapAlleles(gt, lg, i)
  return(gt)
}


#' Find the index of a genotype in a list of genotypes
#'
#' @param gt Genotype (in list format)
#' @param genotypes List of genotypes
#'
#' @keywords internal
#'
#' @return Index of genotype within list, or \code{NA} if the genotype is not found
#'
matchGenotypes <- function(gt, genotypes) {
  m <- NA
  for(i in 1:length(genotypes)) {
    if(identical(gt, genotypes[[i]])) {
      m <- i
      break
    }
  }
  return(m)
}


# A vector containing all possible genotypes in string format
#
# @param genopheno A genopheno object containing all information about a genetic setup.
#
# @return A vector of strings, each representing a possible genotype.
#
getAllGenoStrings <- function(genopheno) {
  # generating table of possible configurations:
  for(lg in 1:length(genopheno$geno)) {
    poss <- list()  # list of possibilities for each locus
    if (genopheno$geno[[lg]]$type == 'autosomal') {
      for(i in 1:genopheno$geno[[lg]]$nloci) {
        locusPoss <- rep(list(1:genopheno$geno[[lg]]$nalleles[i]), 2)
        names(locusPoss) <- c(paste0("LG", lg, "Loc", i, "Pos", 1), paste0("LG", lg, "Loc", i, "Pos", 2))
        poss <- c(poss, locusPoss)
      }
      lgconf <- expand.grid(poss)
    } else if (genopheno$geno[[lg]]$type %in% c('X-linked', 'Z-linked', 'Y-linked', 'W-linked')) {
      nloci <- genopheno$geno[[lg]]$nloci  # for brevity only
      for(i in 1:nloci) {
        locusPoss <- rep(list(0:genopheno$geno[[lg]]$nalleles[i]), 2)
        names(locusPoss) <- c(paste0("LG", lg, "Loc", i, "Pos", 1), paste0("LG", lg, "Loc", i, "Pos", 2))
        poss <- c(poss, locusPoss)
      }
      lgconf <- expand.grid(poss)
      # removing configurations where allele is missing on a chromosome at one but not another locus:
      if (nloci > 1) {
        lgconf <- lgconf[ !someButNotAllColumnZero(lgconf[, 2*(1:nloci) -1]),]
        lgconf <- lgconf[ !someButNotAllColumnZero(lgconf[, 2*(1:nloci)]),]
      }
    } else if (genopheno$geno[[lg]]$type %in% c('maternal', 'paternal')) {
      for(i in 1:genopheno$geno[[lg]]$nloci) {
        locusPoss <- list(1:genopheno$geno[[lg]]$nalleles[i])
        names(locusPoss) <- paste0("LG", lg, "Loc", i)
        poss <- c(poss, locusPoss)
      }
      lgconf <- expand.grid(poss)
    } else {
      stop("Linkage group type not supported.", call. = FALSE)
    }
    if (lg == 1) conf <- lgconf
    else conf <- merge(conf, lgconf)
  }

  genoStrings <- rep("", nrow(conf))
  for(k in 1:nrow(conf)) {
    counter <- 0  # counter for columns in conf dataframe
    for(lg in 1:length(genopheno$geno)) {
      if (genopheno$geno[[lg]]$type %in% c('autosomal', 'X-linked', 'Z-linked', 'Y-linked', 'W-linked')) {
        ploidy <- 2
      } else if (genopheno$geno[[lg]]$type %in% c('maternal', 'paternal')) {
        ploidy <- 1
      } else {
        stop("Linkage group type not supported.", call. = FALSE)
      }
      if(lg>1) genoStrings[k] <- paste0(genoStrings[k], " | ")
      if (genopheno$geno[[lg]]$type %in% c('autosomal', 'maternal', 'paternal')) {
        for(i in 1:genopheno$geno[[lg]]$nloci) {
          if(i>1) genoStrings[k] <- paste0(genoStrings[k], "~")
          for(j in 1:ploidy) {
            counter <- counter + 1
            genoStrings[k] <- paste0(genoStrings[k], genopheno$geno[[lg]]$alleleNames[[i]][ conf[k,counter] ])
          }
        }
      } else if (genopheno$geno[[lg]]$type %in% c('X-linked', 'Z-linked', 'Y-linked', 'W-linked')) {
        for(i in 1:genopheno$geno[[lg]]$nloci) {
          if(i>1) genoStrings[k] <- paste0(genoStrings[k], "~")
          for(j in 1:ploidy) {
            counter <- counter + 1
            genoStrings[k] <- paste0(genoStrings[k], genopheno$geno[[lg]]$alleleNames[[i]][ conf[k,counter] + 1 ])
          }
        }
      } else {
        stop("Linkage group type not supported.", call. = FALSE)
      }
    }
  }
  return(genoStrings)
}


# Determines which genotypes are equivalent to a given genotype.
#
# @param gt a genotype in list format.
# @param genopheno a genopheno object containing all information about a genetic setup.
# @param equivalent Argument determining which genotypes are to be treated as equivalent.
# Possible options are "phase" (genotypes that differ only by phase are equivalent,
# i.e. the order of alleles within a locus does not matter),
# "origin" (phase is important but not the origin of alleles),
# and "none" (no two genotypes are equivalent).
#
# @return A list of genotypes (in list format) that are all equivalent to each other,
#         including the genotype supplied as argument.
#
getEquivalents <- function(gt, genopheno, equivalent) {
  lgs <- length(genopheno$geno)
  if (equivalent == "phase") {
    totalLoci <- sum(sapply(genopheno$geno, function(x) x$nloci))
    swp <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), totalLoci)))
    equis <- vector("list", nrow(swp))
    for(k in 1:nrow(swp)) {
      counter <- 0
      gtEq <- gt
      for(lg in 1:lgs) {
        if (genopheno$geno[[lg]]$type == 'autosomal') {
          for(i in 1:genopheno$geno[[lg]]$nloci) {
            counter <- counter + 1
            if(swp[k, counter]) gtEq <- swapAlleles(gtEq, lg, i)
          }
        } else if (genopheno$geno[[lg]]$type %in% c('X-linked', 'Z-linked', 'Y-linked', 'W-linked')) {
          if ((gt[[lg]][[1]][1] == 0) | (gt[[lg]][[1]][2] == 0)) { # allele missing at first locus (and hence all loci)
            if (all(swp[k, (counter+1):(counter + genopheno$geno[[lg]]$nloci)]))
              gtEq <- swapAllAlleles(gtEq, lg)
            counter <- counter + genopheno$geno[[lg]]$nloci
          } else {
            for(i in 1:genopheno$geno[[lg]]$nloci) {
              counter <- counter + 1
              if(swp[k, counter]) gtEq <- swapAlleles(gtEq, lg, i)
            }
          }
        } else if (genopheno$geno[[lg]]$type %in% c('maternal', 'paternal')) {
            counter <- counter + genopheno$geno[[lg]]$nloci
        } else {
          stop("Linkage group type not recognised.")
        }
      }
      equis[[k]] <- gtEq
    }
  } else if (equivalent == "origin") {
    swp <- expand.grid(rep(list(c(FALSE, TRUE)), lgs))
    equis <- vector("list", nrow(swp))
    for(k in 1:nrow(swp)) {
      gtEq <- gt
      for(lg in 1:lgs) {
        if (genopheno$geno[[lg]]$type %in% c('autosomal', 'X-linked', 'Z-linked', 'Y-linked', 'W-linked')) {
          if(swp[k, lg]) gtEq <- swapAllAlleles(gtEq, lg)
        } else if (genopheno$geno[[lg]]$type %in% c('maternal', 'paternal')) {
          # do nothing
        } else {
          stop("Linkage group type not recognised.")
        }
      }
      equis[[k]] <- gtEq
    }
  } else if (equivalent == "none") {
    equis <- list(gt)
  }
  return(unique(equis))
}

#' Checks if a genotype is a non-canonical genotype
#'
#' With some types of linkage groups such as 'X-linked', there are genotypes that do not
#' usually arise (e.g., absence of both X chromosomes) and should therefore not be displayed by default.
#' This function checks if a genotype is such a non-canonical genotype.
#'
#' @param gt a genotype in list format
#' @param genopheno a genopheno object containing all information about a genetic setup.
#'
#' @return \code{TRUE} or \code{FALSE}
#' @keywords internal
#'
isNoncanonical <- function(gt, genopheno) {
  nc <- FALSE
  for(lg in 1:length(genopheno$geno)) {
    if ((genopheno$geno[[lg]]$type == 'X-linked') & (gt[[lg]][[1]][1] == 0))  # no allele at maternal position
      nc <- TRUE
    if ((genopheno$geno[[lg]]$type == 'Z-linked') & (gt[[lg]][[1]][2] == 0))  # no allele at paternal position
      nc <- TRUE
    if ((genopheno$geno[[lg]]$type == 'Y-linked') & (gt[[lg]][[1]][1] != 0))  # allele at maternal position
      nc <- TRUE
    if ((genopheno$geno[[lg]]$type == 'W-linked') & (gt[[lg]][[1]][2] != 0))  # allele at paternal position
      nc <- TRUE
  }
  return(nc)
}


#' Creates an empty genotype in list format
#'
#' The generated genotype object has the specified number of linkage groups, loci etc.,
#' but only has NAs at all allelic positions.
#'
#' @param genopheno The genetic system that the genotype is part of
#'
#' @keywords internal
#'
#' @return A genotype object in list format
#'
createEmptyGenotype <- function(genopheno) {
  nLGs <- length(genopheno$geno)  # expected number of linkage groups
  geno <- vector("list", nLGs)
  for(lg in 1:nLGs) {  # loop over all linkage groups
    if (genopheno$geno[[lg]]$type %in% c('autosomal', 'X-linked', 'Z-linked', 'Y-linked', 'W-linked')) {
      ploidy <-2
    } else if (genopheno$geno[[lg]]$type %in% c('maternal', 'paternal')) {
        ploidy <- 1
    } else {
      stop("Linkage group not supported")
    }
    geno[[lg]] <- vector("list", genopheno$geno[[lg]]$nloci)
    for(i in 1:genopheno$geno[[lg]]$nloci) { # loop over all loci
      geno[[lg]][[i]] <- rep(NA, ploidy)
    }
  }
  return(geno)
}


# Check if argument is a valid genotype string
#
# @param gtString A genotype in string format.
# @param genopheno genopheno object specifying the genetic setup.
#
# @return Returns TRUE if gtString is a genotype within the genopheno object
#
isValidGenotype <- function(gtString, genopheno) {
  gtString <- tidyGenoString(gtString)
  if (is.null(genopheno$pheno)) {
    if (gtString %in% getAllGenoStrings(genopheno)) return(TRUE)
    else return(FALSE)
  } else {
    if (gtString %in% rownames(genopheno$pheno)) return(TRUE)
    else return(FALSE)
  }
}

#' Helper function to check for which rows of a dataframe some but not all columns are zero
#'
#' @param df A dataframe
#' @keywords internal
#' @return A vector with TRUE and FALSE elements
#'
someButNotAllColumnZero <- function(df) {
  apply(df, 1, function(v) { any(v == 0) & (!all(v == 0)) })
}


