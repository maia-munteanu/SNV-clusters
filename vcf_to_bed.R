library(gtools)

args <- commandArgs(TRUE)
parseArgs <- function(x) {
  res = strsplit(sub("^--", "", x), "=")
  if(length(unlist(res))==1) res[[1]][2]=""
  return(res)
}

argsL <- as.list(do.call("cbind", parseArgs(args))[c(F,T)])
names(argsL) <- as.list(do.call("cbind", parseArgs(args))[c(T,F)])
args <- argsL;rm(argsL)

if(! is.null(args$help)) {
  cat("
      Mandatory arguments:
      --VCF                        - Input VCF to compute the BED file
      --close                      - Maximum distance between the SV breakpoint and the farthest SNV in the low confidence clusters
      --closer                     - Maximum distance between the SV breakpoint and the farthest SNV in the high confidence clusters  
      Optional arguments:
      --output                     - Output folder, otuputs to the current directory by default
      
      --help \n\n")
  q(save="no")
}


if(is.null(args$VCF)) {stop("Option --VCF should be provided")} else{VCF=args$VCF}
if(is.null(args$close)) {stop("Option --close threshold should be provided")} else{close=as.numeric(args$close)}
if(is.null(args$closer)) {stop("Option --closer threshold should be provided")} else{closer=as.numeric(args$closer)}
if(is.null(args$output)) {output="./"} else{output=args$output}

sv_chr = system(paste("bcftools query -f '%CHROM\n' ", VCF, sep=""), intern = T)
sv_pos = as.numeric(system(paste("bcftools query -f '%POS\n' ", VCF, sep=""), intern = T))

closer_df=data.frame(CHR=sv_chr, START=sv_pos-closer, END=sv_pos+closer)
closer_df[which(closer_df$START < 0), "START"] = 1

close_low_df=data.frame(CHR=sv_chr, START=sv_pos-close, END=sv_pos-closer)
close_low_df=close_low_df[which(close_low_df$END>0),]
close_high_df=data.frame(CHR=sv_chr, START=sv_pos+closer, END=sv_pos+close)
close_df=rbind(close_low_df,close_high_df)

chrOrder <-c((1:22),"X","Y")
close_df$CHR <- factor(close_df$CHR, chrOrder)
close_df=close_df[order(close_df$CHR,close_df$START),]

cluster_df=data.frame(CHR=sv_chr, START=sv_pos-close, END=sv_pos+close)
cluster_df[which(cluster_df$START<0),"START"]=1


sample_name=strsplit(VCF,"[.]")[[1]][1]
cat("\n \n","Outputting beds for sample ",sample_name,"\n")

setwd(output)
write.table(closer_df,file=paste(sample_name,"_0_",closer/1000,"kb_closer.bed",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
write.table(close_df,file=paste(sample_name,"_",closer/1000,"kb_",close/1000,"kb_close.bed",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
write.table(cluster_df, file=paste(sample_name,"_0_",close/1000,"kb_cluster.bed",sep=""), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")


