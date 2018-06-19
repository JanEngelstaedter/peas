#' Generates a table of all possible egg configurations
#'
#' In this table, each egg configuration is a way in which alleles could be drawn from a
#' genotype to produce an egg haplotype. The table also has a column with expected
#' fractions of these configurations.
#' Given an actual genotype, each configuration can then be used to generate the actual egg haplotype.
#'
#' @param genopheno A genopheno object containing all information about the genetic setup.
#'
#' @return A data frame whose first columns represent all loci.
#' For each locus, a number 1 represents the maternal allele, 2 represents the paternal allele,
#' and 0 represents no allele. The last column contains the expected fraction of each configuration.
#'
getEggConfigurations <- function(genopheno) {
  # generating table of possible configurations:
  poss <- list()  # list of possibilities for each locus
  for(lg in 1:length(genopheno$geno)) {
    if (genopheno$geno[[lg]]$type %in% c('autosomal', 'XY')) {
      poss <- c(poss, rep(list(c(1,2)), genopheno$geno[[lg]]$nloci))
    } else if (genopheno$geno[[lg]]$type = 'ZW') {
      poss <- c(poss, rep(list(c(1,0)), genopheno$geno[[lg]]$nloci))
    } else if (genopheno$geno[[lg]]$type = 'maternal') {
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

#' Generates a table of all possible sperm configurations
#'
#' In this table, each sperm configuration is a way in which alleles could be drawn from a
#' genotype to produce an egg haplotype. The table also has a column with expected
#' fractions of these configurations.
#' Given an actual genotype, each configuration can then be used to generate the actual sperm haplotype.
#'
#' @param genopheno A genopheno object containing all information about the genetic setup.
#'
#' @return A data frame whose first columns represent all loci.
#' For each locus, a number 1 represents the maternal allele, 2 represents the paternal allele,
#' and 0 represents no allele. The last column contains the expected fraction of each configuration.
#'
getSpermConfigurations <- function(genopheno) {
  # generating table of possible configurations:
  poss <- list()  # list of possibilities for each locus
  for(lg in 1:length(genopheno$geno)) {
    if (genopheno$geno[[lg]]$type %in% c('autosomal', 'ZW')) {
      poss <- c(poss, rep(list(c(1,2)), genopheno$geno[[lg]]$nloci))
    } else if (genopheno$geno[[lg]]$type = 'XY') {
      poss <- c(poss, rep(list(c(1,0)), genopheno$geno[[lg]]$nloci))
    } else if (genopheno$geno[[lg]]$type = 'maternal') {
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
