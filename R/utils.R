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
#'
#' @examples
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
  for(p in pos) substr(strings,p,p)<-"_"
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
#' @return A genotype object in list format
#'
convertGenoStringToList <- function(genoString, genopheno) {
  genoString <- tidyGenoString(genoString)  # removing excess spaces etc.
  nLGs <- length(genopheno$geno)  # expected number of linkage groups
  genoList <- vector("list", nLGs)
  spos <- 1  # position within genoString
  for(lg in 1:nLGs) {  # loop over all linkage groups
    if (genopheno$geno[[lg]]$type == 'autosomal') {
      ploidy <-2
    } else {
      stop("Non-autosomal linkage groups not yet implemented")
    }
    if (lg > 1) spos <- spos + 3 # skipping the " | " seperator
    genoList[[lg]] <- vector("list", genopheno$geno[[lg]]$nloci)
    for(i in 1:genopheno$geno[[lg]]$nloci) { # loop over all loci
      len <- nchar(genopheno$geno[[lg]]$alleleNames[[i]][1])
      if (i>1) spos <- spos + 1  # skip tilde
      genoList[[lg]][[i]] <- rep(NA, ploidy)
      for(pos in 1:ploidy) { # loop over the alleles per locus
        for(j in 1:genopheno$geno[[lg]]$nalleles[i]) { # loop over allele names
          if (substr(genoString, spos, spos + len -1) == genopheno$geno[[lg]]$alleleNames[[i]][j]) {  # match!
            genoList[[lg]][[i]][pos] <- j
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
#' @return Genotype in string format
#'
#' @examples
convertGenoListToString <- function(genoList, genopheno) {
  genoString <- "" # start with empty string
  for(lg in 1:length(genoList)) { # loop over all linkage groups
    if(lg>1) genoString <- paste0(genoString, " | ")
    for(i in 1:length(genoList[[lg]])) { # loop over all loci
      if(i>1) genoString <- paste0(genoString, "~")
      for(j in 1:length(genoList[[lg]][[i]]))  # loop over alleles at locus
        genoString <- paste0(genoString, genopheno$geno[[lg]]$alleleNames[[i]][[ genoList[[lg]][[i]][j] ]])
    }
  }
  return(genoString)
}

#' Swapping alleles at the two positions of a locus
#'
#' @param gt Genotype in list format
#' @param lg Linkage group numer
#' @param loc Locus number
#'
#' @return Genotype in list format
#'
#' @examples
swapAlleles <- function(gt, lg, loc) {
  dummy <- gt[[lg]][[loc]][1]
  gt[[lg]][[loc]][1] <- gt[[lg]][[loc]][2]
  gt[[lg]][[loc]][2] <- dummy
  return(gt)
}


#' Swapping the two alleles at all loci of a linkage group
#'
#' @param gt Genotype in list format.
#' @param lg Linkage group where alleles should be swapped.
#'
#' @return A genotype in list format.
#'         This genotype is identical to the supplied one except for the swapped linkage group.
#'
#' @examples
swapAllAlleles <- function(gt, lg) {
  for(i in 1:length(gt[[lg]])) gt <- swapAlleles(gt, lg, i)
  return(gt)
}


#' Find the index of a genotype in a list of genotypes
#'
#' @param gt Genotype (in list format)
#' @param genotypes List of genotypes
#'
#' @return Index of genotype within list, or \code{NA} if the genotype is not found
#'
#' @examples
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


#' A vector containing all possible genotypes in string format
#'
#' @param genopheno A genopheno object containing all information about a genetic setup.
#'
#' @return A vector of strings, each representing a possible genotype.
#'
getAllGenoStrings <- function(genopheno) {
  # generating table of possible configurations:
  poss <- list()  # list of possibilities for each locus
  for(lg in 1:length(genopheno$geno)) {
    if (genopheno$geno[[lg]]$type == 'autosomal') {
      for(i in 1:genopheno$geno[[lg]]$nloci)
        poss <- c(poss, rep(list(1:genopheno$geno[[lg]]$nalleles[i]), 2))
    } else {
      stop("Non-autosomal linkage groups not yet implemented.")
    }
  }
  conf <- expand.grid(poss)

  genoStrings <- rep("", nrow(conf))
  for(k in 1:nrow(conf)) {
    counter <- 0  # counter for columns in conf dataframe
    for(lg in 1:length(genopheno$geno)) {
      if(lg>1) genoStrings[k] <- paste0(genoStrings[k], " | ")
      if (genopheno$geno[[lg]]$type == 'autosomal') {
        for(i in 1:genopheno$geno[[lg]]$nloci) {
          if(i>1) genoStrings[k] <- paste0(genoStrings[k], "~")
          for(j in 1:2) {
            counter <- counter + 1
            genoStrings[k] <- paste0(genoStrings[k], genopheno$geno[[lg]]$alleleNames[[i]][ conf[k,counter] ])
          }
        }
      } else {
        stop("Non-autosomal linkage groups not yet implemented.")
      }
    }
  }
  return(genoStrings)
}


#' Determines which genotypes are equivalent to a given genotype.
#'
#' @param gt A genotype in list format.
#' @param genopheno A genopheno object containing all information about a genetic setup.
#' @param equivalent Argument determining which genotypes are to be treated as equivalent.
#' Possible options are "phase" (genotypes that differ only by phase are equivalent,
#' i.e. the order of alleles within a locus does not matter),
#' "origin" (phase is important but not the origin of alleles),
#' and "none" (no two genotypes are equivalent).
#'
#' @return A list of genotypes (in list format) that are all equivalent to each other,
#'         including the genotype supplied as argument.
#'
getEquivalents <- function(gt, genopheno, equivalent) {
  lgs <- length(genopheno$geno)
  if (equivalent == "phase") {
    totalLoci <- sum(sapply(genopheno$geno, function(x) x$nloci))
    swp <- expand.grid(rep(list(c(FALSE, TRUE)), totalLoci))
    equis <- vector("list", nrow(swp))
    for(k in 1:nrow(swp)) {
      counter <- 0
      gtEq <- gt
      for(lg in 1:lgs) {
        if (genopheno$geno[[lg]]$type %in% c('autosomal', 'XY', 'ZW')) {
          for(i in 1:genopheno$geno[[lg]]$nloci) {
            counter <- counter + 1
            if(swp[k, counter]) gtEq <- swapAlleles(gtEq, lg, i)
          }
        } else if (genopheno$geno[[lg]]$type == 'maternal') {
          # do nothing
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
        if (genopheno$geno[[lg]]$type %in% c('autosomal', 'XY', 'ZW')) {
          if(swp[k, lgs]) gtEq <- swapAllAlleles(gtEq, lg)
        } else if (genopheno$geno[[lg]]$type == 'maternal') {
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



# The functions below should no longer be necessary:

# swap two characters in a string s at positions pos1 and pos2:
#swap.chars<-function(s,pos1,pos2)
#{
#  c<-substr(s,pos1,pos1)
#  substr(s,pos1,pos1)<-substr(s,pos2,pos2)
#  substr(s,pos2,pos2)<-c
#  return(s)
#}

# insert placeholders ("_") in all genotype names at positions
# specified by the vector pos

#insert.placeholders<-function(pos)
#{
#  genotypes<-names(PHEN)
#  for(p in pos) substr(genotypes,p,p)<-"_"
#  return(genotypes)
#}


