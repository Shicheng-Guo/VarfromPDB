extract_clinvar <-
function(keyword, localPDB.path=paste(getwd(), "localPDB",sep="/"), type="both",
        HPO.disease=NULL, genelist=NULL){
    morbidmap=paste(localPDB.path,"morbidmap.txt",sep="/")
    if(file.exists(localPDB.path)){
         if(file.exists(paste(localPDB.path,"variant_summary.txt.gz",sep="/"))){
             clinvar <- paste(localPDB.path,"variant_summary.txt.gz",sep="/")
             }else{
                 clinvar <- NULL
         }        

         if(file.exists(paste(localPDB.path,"gene_condition_source_id",sep="/"))){
             gene2dis <- paste(localPDB.path,"gene_condition_source_id",sep="/")
             }else{
                 gene2dis <- NULL
         }        

         if(file.exists(paste(localPDB.path,"NAMES.csv.gz",sep="/"))){
             medgene.names <- paste(localPDB.path,"NAMES.csv.gz",sep="/")
             }else{
                 medgene.names <- NULL
         }        

         if(file.exists(paste(localPDB.path,"GRtitle_shortname_NBKid.txt",sep="/"))){
             genereview <- paste(localPDB.path,"GRtitle_shortname_NBKid.txt",sep="/")
             }else{
                 genereview <- NULL
         }        

        }else{
             clinvar <- NULL
             gene2dis <- NULL
             medgene.names <- NULL
             genereview <- NULL
    }     

    #check clinvar database
    if(is.null(clinvar)){
       clinvar <- "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
       download.path <- paste(getwd(),"localPDB",sep="/")
       if(!file.exists(download.path))
          dir.create(download.path )
       options(timeout = 300)
       if( !file.exists(paste(download.path,"variant_summary.txt.gz",sep="/")))
           download.file(clinvar,
               paste(download.path,"variant_summary.txt.gz",sep="/"),method="auto")
       clinvar <- paste(download.path,"variant_summary.txt.gz",sep="/")
    }

    if(substr(clinvar,nchar(clinvar)-1,nchar(clinvar)) == "gz"){
        clinvar <- read.delim(gzfile(clinvar))
        }else{
            clinvar <- read.delim(clinvar)
    }       

    if(is.null(gene2dis)){
       gene2dis <- "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id"
       download.path <- paste(getwd(),"localPDB",sep="/")
       if(!file.exists(download.path))
           dir.create(download.path )
       options(timeout = 300)
       if(!file.exists(paste(download.path,"gene_condition_source_id",sep="/")) )
           download.file(gene2dis,paste(download.path,"gene_condition_source_id",sep="/"),method="auto")
       gene2dis <- paste(download.path,"gene_condition_source_id",sep="/")
    }
       gene2dis <- read.delim(gene2dis)
       

   ## HPO
       if(is.null(HPO.disease)){
          HPO.disease.check <- pheno_extract_HPO(keyword= keyword)
          HPO.disease <- as.character(unique(HPO.disease.check[grep("OMIM",HPO.disease.check[,1]),1]))
       }


   ##function: extract the phenotypeID from the colnames of PhenotypeIDS in clinvar_summary file
   ##databaseID="MedGen","OMIM","GeneReviews",unavailable: "SNOMED CT","Orphanet",
   extract.PhenotypeID <- function(x,database){
      #x=clinvar.extr[49,"PhenotypeIDS"];database="OMIM"
       x.split <- setdiff(unlist(strsplit(as.character(x),";")),"")
       ids <- grep(database,x.split,ignore.case = TRUE)
       if(length(ids)==0){
          id.extr <- ""
          }else if(length(ids)==1){
             x.split.2 <- unlist(strsplit(x.split[ids],","))
             id.extr <- x.split.2[grep(database,x.split.2,ignore.case = TRUE)]
             if(length(id.extr)>1){
                  id.extr <- unique(id.extr)
                  if(length(id.extr)>1){
                      id.extr <- paste(id.extr,collapse=",")
                      }
                  }
             }else{
              id.extr <- c()
              for(j in ids){
                  x.split.2 <- unlist(strsplit(x.split[j],","))
                  id.extr.j <- x.split.2[grep(database,x.split.2,ignore.case = TRUE)]
                  id.extr <- c(id.extr,id.extr.j)
              }  
              id.extr <- paste(id.extr,collapse=",")
       }     
       return(id.extr) 
   }
##########################

 ##begin to search 
       if(!is.null(keyword)){      
          gene2dis.d <- gene2dis[grep_split(keyword,gene2dis[,"DiseaseName"]),]
          pheno.yes <- as.character(gene2dis.d[,"DiseaseName"])
          }else if((is.null(keyword)& !is.null(HPO.disease)) | (is.null(keyword) & !is.null(genelist))){
            gene2dis.d <- c()
            pheno.yes <- c()  
            }else{
              gene2dis.d <- gene2dis
              pheno.yes <- c()
       }           
     
       if(!is.null(HPO.disease)){
          HPO.disease.no <- unlist(lapply(HPO.disease,function(x) unlist(strsplit(x,"OMIM:"))[2]))
          gene2dis.d2 <- gene2dis[is.element(gene2dis[,"DiseaseMIM"],HPO.disease.no),]
          pheno.yes2 <- as.character(gene2dis.d2[,"DiseaseName"])
          gene2dis.d <- rbind(gene2dis.d,gene2dis.d2)
          pheno.yes <- union(pheno.yes,pheno.yes2)
       }
   
   #for a given genelist
       if(!is.null(genelist)){
          gene2dis.d3 <- gene2dis[is.element(gene2dis[,"GeneSymbol"],genelist),]
          gene2dis.d <- rbind(gene2dis.d,gene2dis.d3)
       }

       gene2dis.extr <- unique(gene2dis.d)   
      gene2dis.extr$pheno.check <- "no"
      gene2dis.extr[is.element(gene2dis.extr$DiseaseName,pheno.yes),"pheno.check"] <-  "yes"
           
      genes <- unique(as.character(gene2dis.extr[,2]))
     mim.id.pheno <- unique(gene2dis.extr[gene2dis.extr[,"pheno.check"] == "yes",7])
     mim.id.pheno <- mim.id.pheno[!is.na(mim.id.pheno)]
     
#      morbidmap <- read.table(morbidmap,header= FALSE,sep="|",quote = "")           
       morbidmap <- read.delim(morbidmap,comment.char = "#")           
      colnames(morbidmap) <- c("disease","gene","gene.mim.no","location")
     
     ##extract the variants in the genes
     clinvar.extr <- clinvar[is.element(clinvar[,"GeneSymbol"],genes),]
     
## extract the variants from summary file directly
      clinvar.d <- clinvar[grep_split(keyword,clinvar[,"PhenotypeList"]),]
     
     clinvar.extr <- unique(rbind(clinvar.extr,clinvar.d))
     extract <- list(gene2dis.extr,clinvar.extr)
     names(extract) <- c("gene2dis","variants")
     return(extract)
}
