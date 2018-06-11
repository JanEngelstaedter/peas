


#' Generates a new genetic system
#'
#' This function generates a new genetic system that stores information about the number of loci, alleles at each locus
#' as well as recombination rates and mode of inheritance.
#'
#' @param nloci The number of loci. The default value is 1 but any (reasonable) number can be chosen.
#' @param nalleles The number of alleles for each locus. This should be a vector of length \code{nloci} containing positive integer numbers.
#' By default, all loci are set to be biallelic.
#' @param alleleNames The names for the alleles at the different loci.
#' The default value is 'standard', which automatically produces allele names 'a' and 'A' for the first locus, 'b' and 'B' for the second locus etc.
#' For other allele names, \code{alleleNames} should be a list of nloci elements, each of which needs to be a vector containing the allele names.
#' Allele names need to be strings (class \code{character}) and can contain letters, numbers, or symbols (except '~', '.' and '|').
#' @param rec A vector of length \code{nloci-1} specifying recombination rates between the loci.
#' @param type This parameter specifies the ploidy and inheritance for the loci.
#' The default value is "autosomal", corresponding to diploid, autosomal inheritance.
#' Other types that are currently supported include 'XY', 'ZW', and 'cytoplasmic'.
#'
#' @return An object of class \code{genopheno} that stores all of the specifications provided about the genetic setup
#' and that can then be extended to also include other linkage groups and corresponding phenotypes.
#' @export
#'
#' @examples
new_genopheno <- function(nloci = 1, nalleles = rep(2, nloci), alleleNames = 'standard', rec = NULL, type = 'autosomal') {

  # initial checking if all arguments make sense and are consistent:

  if ((!is.numeric(nloci)) | (length(nloci)>1) | (nloci < 1))
    stop("Number of loci must be a positive integer number.")
  if (length(nalleles) != nloci)
    stop(paste0("Vector of length nloci=",nloci," expected for argument 'nalleles'."))
  if (alleleNames[1] != 'standard') {
    if ((!is.list(alleleNames)) | (length(alleleNames)!=nloci))
      stop(paste0("List of length nloci=", nloci, " expected for argument 'alleleNames'."))
    if (!all.equal(sapply(alleleNames, length), nalleles))
      stop("Number of alleles does not correspond to allele names.")
  }
  if ((nloci == 1) & (!is.null(rec)))
    stop("No recombination rates can be set when there is only one locus.")
  if (length(rec) != (nloci-1))
    stop(paste("Vector of lenght",nloci-1,"expected for argument 'rec'."))
  if (!(type %in% c('autosomal')))
    stop("Inheritance specified by 'type' not supported.")

  if (alleleNames[1]=='standard')
  {
    alleleNames<-vector("list",NLOCI)
    if (max(NALLELES)==2)
      for (locus in 1:NLOCI) alleleNames[[locus]]<-c(letters[locus],LETTERS[locus])
      else
        for (locus in 1:NLOCI) alleleNames[[locus]]<-1:NALLELES[locus]
  }

  gp <- list(geno = list(LG1 = list(nloci = as.integer(nloci),
                              nalleles = as.integer(nalleles),
                              alleleNames = alleleNames,
                              rec = rec,
                              type = type)),
       pheno = NULL)
  class(gp) <- "genopheno"
  gp
}


allelicStates<-list()
for (locus in 1:NLOCI)
{
  for (pos in 1:2)
    allelicStates[[3*(locus-1)+pos]]<-alleleNames[[locus]]
  if (locus<NLOCI)
    allelicStates[[3*locus]]<-"~"
}
genotypes<-sort(apply(expand.grid(allelicStates),1,paste,collapse=''))
PHEN<<-rep(NA,length(genotypes))
names(PHEN)<<-genotypes
