
options(stringsAsFactors = FALSE)
setwd("/g/strcombio/fsupek_cancer1/SV_clusters_project/")

#read metadata file
metadata=read.table(file = "metadata.tsv", sep = '\t', header = TRUE)
#find all cancer types in the datset
cancer_types=unique(metadata[,"primaryTumorLocation"])
#print cancer types
cat("\n","\n","#################################### \n","############### INFO ############### \n","\n")
cat("The dataset contains ",length(cancer_types), " different cancer types \n")
cat(cancer_types,sep="\n")
cat("\n","\n","\n","##################################### \n","############### INPUT ############### \n","\n")

#Ask for input - no. of different cancer types and the names of the cancer types
cat("Enter the number of different cancer types to extract: \n")
no_of_cancer_types=as.numeric(readLines(con="stdin", 1))
cat("Enter ", no_of_cancer_types, " cancer type name/s below: \n")

#Ask for input - all cancer types requested
query_vector=vector(mode="character", length=no_of_cancer_types)
for (query_no in 1:no_of_cancer_types){
  incorrect=0
  while(incorrect==0){ 
    temp=readLines(con="stdin",1)
    if (temp %in% cancer_types){
       if (temp %in% query_vector){cat("Cancer type has been entered before, try a different cancer type: \n")}
       else{
           query_vector[query_no]=temp
           incorrect=1}}      
    else{cat("Cancer type entered does not appear in this dataset, try again: \n")}
}}

cat("\n","\n")

#list with samples for each cancer type entered
query_list=vector(mode = "list", length = no_of_cancer_types)
for (query_no in 1:no_of_cancer_types){
    query_list[[query_no]]=metadata[which(metadata$primaryTumorLocation==query_vector[query_no]),1]}
names(query_list)=query_vector
cat("List of samples from queried cancer type/s found in the metadata file:","\n")
print(query_list)

cat("\n","\n","###################################### \n","############### OUTPUT  ############## \n","\n")
cat("Searching the Hartwig directory on agendas for chosen samples","\n")

#create a dataframe with all the samples currently available on Agendas
hartwig_files = list.files("/g/strcombio/fsupek_cancer1/HARTWIG_vcfs/somatic_update_Nov2021", full.names = T)
hartwig_files = hartwig_files[which(!grepl("tsv", hartwig_files))]
hartwig_files = hartwig_files[which(!grepl("json", hartwig_files))]
for(file in hartwig_files){  
  # get the files
  sample_name=strsplit(file,"/")[[1]][7]
  snvs = paste(file, "/purple/", sample_name, ".purple.somatic.vcf.gz", sep="")
  svs = paste(file, "/purple/", sample_name, ".purple.sv.vcf.gz", sep="")
  if (file.exists(snvs)==TRUE && file.exists(svs)==TRUE){
     temp_df = data.frame(sample = sample_name, sv = svs, snv = snvs) }
  if(match(file, hartwig_files) == 1) {df = temp_df} else {df = rbind(df, temp_df)}
}

#Match the samples from the metadata file for each cancer type with the list of samples currently available and create new database 
for (i in 1:no_of_cancer_types){
    temp_df=df[match(query_list[[i]], df$sample),]
    temp_df=temp_df[which(!is.na(temp_df$sample)),]
    cat(length(temp_df[,1]), "/", length(query_list[[i]])," ",query_vector[i], " samples found in the given directory from total of samples in the metadata file \n")
    if (i==1){merge_df=temp_df} else {merge_df=rbind(merge_df,temp_df)}
}


#write the .txt file 
merge_names=paste(query_vector,collapse="+")
merge_names=gsub(" ","-",merge_names)
merge_names=gsub("/","-",merge_names)
write.table(merge_df, file=paste("Hartwig_",merge_names,"_samples.csv",sep=""), quote = F, row.names = F, col.names = T, sep = "\t")

cat("\n","Returned a .csv file with  ",length(merge_df[,1]),"  samples","\n","\n","\n")

