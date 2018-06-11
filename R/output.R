

#' Print information about genotypes and corresponding phenotypes.
#'
#' @param gp Object of class \code{genopheno}
#' @param what What information should be printed.
#' With the default value 'all', both information about the genetic system and the phenotypes is printed.
#' With values 'geno' and 'pheno', only information about the genotypes or phenotypes is printed, respectively.
#'
#' @return NULL
#' @export
#'
#' @examples
print.genopheno <- function(gp, what = 'all') {

  if (!(what %in% c('all', 'geno', 'pheno')))
    stop("Invalid value for argument 'what'.")

  if ((what == 'all') | (what == 'geno')) {
    nLGs <- length(gp$geno) # number of linkage groups
    if (nLGs == 1)
      cat("Genetic system comprising one linkage group\n")
    else
      cat("Genetic system comprising",nLGs, "linkage groups\n", sep="")
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
  }
#  if ((what == 'all') | (what == 'geno')) {
#
#  }
}
