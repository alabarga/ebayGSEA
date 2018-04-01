#### ebayGSEA.R

#### Description: implements a GSEA designed for Illumina Infinium Methylation beadchips, based on an empirical Bayes method to rank genes based on their level of differential methylation, subsequently assessing enrichment of biological terms using this ranked list.
#### Author: Andrew Teschendorff (a.teschendorff@ucl.ac.uk)
#### Date: 1st Apr.2018

### load libraries
library(parallel);
library(globaltest);
library(org.Hs.eg.db);
library(kpmt);

### Auxiliary functions
convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )

  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }

  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

PasteVector <- function(v){

  vt <- v[1];
  if(length(v) > 1){
   for(g in 2:length(v)){
    vt <- paste(vt,v[g],sep=" ")

   }
  }
  vt <- paste(vt," EnD",sep="");
  out.v <- sub(" EnD","",vt);
  out.v <- sub("NA , ","",out.v);
  out.v <- sub(" , NA","",out.v);
  out.v <- sub(" , NA , "," , ",out.v);
  return(out.v);
}


#### Function to rank genes on Illumina platform using the empirical Bayes global test 
doGT <- function(pheno.v,data.m,model=c("linear"),array=c("450k","850k")){

#### load objects
load("dualmap450kEID.Rd");
load("dualmap850kEID.Rd");

if(array=="450k"){
   subsets <- lapply(mapEIDto450k.lv,intersect,rownames(data.m));
}
else {
   subsets <- lapply(mapEIDto850k.lv,intersect,rownames(data.m));
}
nrep.v <- unlist(lapply(subsets,length));
selG.idx <- which(nrep.v>0);    
gt.o <- gt(response=pheno.v,alternative=t(data.m),model=model,directional = FALSE, standardize = FALSE, permutations = 0, subsets=subsets[selG.idx],trace=TRUE);
resGT.m <- as.matrix(result(gt.o));
tmp.s <- sort(resGT.m[,2],decreasing=TRUE,index.return=TRUE);
sresGT.m <- resGT.m[tmp.s$ix,];
return(sresGT.m);
}


### Functions to perform GSEA using wilcox-test and the known population median test (threshold independent)
gseaWTfn <- function(termEID.v,rankEID.v,minN=5){

    commonEID.v <- intersect(termEID.v,rankEID.v);
    nrep <- length(commonEID.v);
    if(length(commonEID.v)>=minN){
    otherEID.v <- setdiff(rankEID.v,termEID.v);
    match(commonEID.v,rankEID.v) -> rank1.idx;
    match(otherEID.v,rankEID.v) -> rank2.idx;
    wt.o <- wilcox.test(rank1.idx,rank2.idx,alt="less");
    pv <- wt.o$p.value;
    n1 <- length(rank1.idx);
    n2 <- length(rank2.idx);
    auc <- 1 - wt.o$stat/(n1*n2);
    ### now do kpmt
    pop.v <- 1:length(rankEID.v);
    names(pop.v) <- rankEID.v;
    obs.v <- commonEID.v;
    pvKPMT <- kpmt(pop=pop.v,obs=obs.v,tail="lower")[[4]];
    out.v <- c(nrep,auc,pv,pvKPMT);
    }
    else {
      out.v <- c(nrep,0,1,1);
    }
    return(out.v);
}

doGSEAwt <- function(rankEID.v,listEZ.lv,ncores=4,minN=5,adjPVth=0.05){

gseaWT.m <- matrix(unlist(mclapply(listEZ.lv,gseaWTfn,rankEID.v,mc.cores=ncores,minN=minN)),ncol=4,byrow=TRUE)
colnames(gseaWT.m) <- c("nREP","AUC","P(WT)","P(KPMT)");
rownames(gseaWT.m) <- names(listEZ.lv);

tmp.s <- sort(gseaWT.m[,3],decreasing=FALSE,index.return=TRUE);
sgseaWT.m <- gseaWT.m[tmp.s$ix,];
padj.v <- p.adjust(sgseaWT.m[,3],method="BH");

sel.idx <- which(padj.v <= adjPVth);
topGSEAwt.lm <- list();
if(length(sel.idx)>1){

topGSEAwt.m <- cbind(sgseaWT.m[sel.idx,],padj.v[sel.idx]);
colnames(topGSEAwt.m) <- c("nREP","AUC","P(WT)","P(KPMT)","adjP");

topGSEAwt.lm[[1]] <- topGSEAwt.m;
tmp.s <- sort(topGSEAwt.m[,2],decreasing=TRUE,index.return=TRUE);
topGSEAwt.lm[[2]] <- topGSEAwt.m[tmp.s$ix,];
names(topGSEAwt.lm) <- c("Rank(P)","Rank(AUC)");
}

else if (length(sel.idx)==1) {
    topGSEAwt.v <- as.vector(c(sgseaWT.m[sel.idx,],padj.v[sel.idx]));
    names(topGSEAwt.v) <- c("nREP","AUC","P(WT)","P(KPMT)","adjP");
    topGSEAwt.lm <- list("Rank(P)"=topGSEAwt.v,"Rank(AUC)"=topGSEAwt.v,"POI"=rownames(sgseaWT.m)[sel.idx]);
}

return(topGSEAwt.lm);
}


### Example
### Assume linear phenotype is in pheno.v, and data matrix in data.m, and that array is Illumina 450k, biological terms annotated to EID are in listEZ.lv;

sgt.m <-doGT(pheno.v,data.m,array=c("450k"));
topGSEA.lm <- doGSEAwt(rankEID.v=rownames(sgt.m),listEZ.lv,ncores=4,minN=5,adjPVth=0.3);

