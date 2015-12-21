\name{genes_compile}
\alias{genes_compile}
\title{
Compile the disease-related genes from multiple public databases
}
\description{
To compile a gene set related to a disease especially for a rare disease from 
multiple databases, including HPO, orphanet, omim, clinvar and uniprot.
}
\usage{
genes_compile(HPO, orphanet, omim, clinvar, uniprot,
       localPDB = paste(getwd(),"localPDB",sep="/"))
}
\arguments{
  \item{HPO}{
the object from \emph{pheno.extract.HPO} function.
}
  \item{orphanet}{
the object from \emph{extract.genes.orphanet} function.
}
  \item{omim}{
the object from \emph{extract.omim} function.
}
  \item{clinvar}{
the object from \emph{extract.clinvar} function.
}
  \item{uniprot}{
the object from \emph{extract.uniprot} function.
}
  \item{localPDB}{
the path of localized public data bases.
}
}
\details{
The relation between genes and a phenotype in different databases can be 
intergrated automatically. 
}
\value{
A matrix containing the following information
   \item{GeneSymbol }{gene symbol in HGNC.}
   \item{chr }{chorosome of a gene.}
   \item{strand }{strand of a gene.}
   \item{start }{start position (hg19) of gene.}
   \item{end }{end position (hg19) of a gene.}
   \item{EntrezGeneID }{Entrez GeneID}
   \item{ApprovedName }{Approved gene name in HGNC.}
   \item{Synonyms }{gene Synonyms.}
   \item{HPO }{the phenotype in HPO.}
   \item{Orphanet }{the phenotype in orphanet.}
   \item{OMIM }{the phenotype in OMIM.}
   \item{ClinVar }{the phenotype in ClinVar.}
   \item{Uniprot }{the phenotype in Uniprot.}
}
\author{
Zongfu Cao
}
\seealso{
 \code{\link{pheno_extract_HPO}}, \code{\link{extract_omim}}, 
 \code{\link{extract_genes_orphanet}}, \code{\link{extract_clinvar}}, 
 \code{\link{extract_uniprot}}
}
\examples{
## compile the gene-disease relationship from multiple databases
#genesPDB <- genes_compile(HPO = HPO.Joubert, orphanet = orphanet.joubert,
#                          omim = genes.omim,
#                          clinvar = genes.clinvar,
#                          uniprot = genes.uniprot)
}                          

\keyword{ phenotype }
\keyword{ gene }
