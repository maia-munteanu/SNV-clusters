#create a file that holds the paths to all the Hartwig samples
setwd("/g/strcombio/fsupek_cancer1/SV_clusters_project/")

hartwig_files = list.files("/g/strcombio/fsupek_cancer1/HARTWIG_vcfs/somatic_update_Nov2021", full.names = T)
hartwig_files = hartwig_files[which(!grepl("tsv", hartwig_files))]
hartwig_files = hartwig_files[which(!grepl("json", hartwig_files))]

for(file in hartwig_files){
  print(paste(date(), " INFO: working on sample number", match(file, hartwig_files)))
  
  # get the files
  sample_name=strsplit(file,"/")[[1]][7]
  snvs = paste(file, "/purple/", sample_name, ".purple.somatic.vcf.gz", sep="")
  svs = paste(file, "/purple/", sample_name, ".purple.sv.vcf.gz", sep="")
  
  if (file.exists(snvs)==TRUE && file.exists(svs)==TRUE){
     temp_df = data.frame(sample = sample_name, sv = svs, snv = snvs) }
  
  if(match(file, hartwig_files) == 1) {df = temp_df} else {df = rbind(df, temp_df)}
}

write.table(df, file="Hartwig_all_samples.csv", quote = F, row.names = F, col.names = T, sep = "\t")
print(paste("Returned a .csv file with ", length(df$sample), " samples"))
