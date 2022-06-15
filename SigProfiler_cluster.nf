#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

// Copyright (C) 2022 IRB Barcelona

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

log.info ""
log.info "-------------------------------------------------------------------------"
log.info "  SNV-clusters: nextflow pipeline to extract SNVs found near SVs "
log.info "          count them and discover mutational signatures   "
log.info "-------------------------------------------------------------------------"
log.info "Copyright (C) IRB Barcelona"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "-------------------------------------------------------------------------"
log.info ""

params.help = null
if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info '  USAGE              '
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run main.nf --input_file Hartwig_all_samples.csv --closer_value 1000 --close_value 10000 --output_folder /g/*/*cancer1*/SV* --fasta_file /g/*/*cancer1*/SV*/hg19.fasta --hg19 /g/*/*cancer1*/SV*/hg19.genome'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --input_file                   FILE           Input .csv file containing 3 columns: sample name, SV vcf path and SNV vcf path.'
    log.info '    --closer_value                 INTEGER        Maximum distance (in bp) between the SV breakpoint and the farthest SNV in the high confidence clusters.'
    log.info '    --close_value                  INTEGER        Maximum distance (in bp) between the SV breakpoint and the farthest SNV in the low confidence clusters.'
    log.info '    --output_folder                FOLDER         Output folder.'
    log.info '    --fasta_ref                    FILE           Fasta reference file.'
    log.info '    --hg19                         FILE           Tab delimited .genome file containing chromosome lengths for hg19.'
    log.info ''
    log.info 'Flags:'
    log.info '    --help                                        Display this message'
    log.info ''
    exit 0
}

params.input_file = "/g/strcombio/fsupek_cancer1/SV_clusters_project/Hartwig_all_samples_filt.csv"
params.closer_value = 2000
params.close_value = 10000
params.output_folder = "/g/strcombio/fsupek_cancer1/SV_clusters_project/SigProfiler"
params.fasta_ref = "/g/strcombio/fsupek_cancer1/SV_clusters_project/hg19.fasta"
params.hg19 = "/g/strcombio/fsupek_cancer1/SV_clusters_project/hg19.genome"


pairs_list = Channel.fromPath(params.input_file, checkIfExists: true).splitCsv(header: true, sep: '\t', strip: true)
                   .map{ row -> [ row.sample, file(row.sv), file(row.snv) ] }.view()
hg19 = file(params.hg19)
fasta_ref=file(params.fasta_ref)

process make_sv_beds {

       publishDir params.output_folder+"/SV_BEDs/", mode: 'copy', pattern: '*.bed'

       input:
       set val(sample), file(sv), file(snv) from pairs_list
       file hg19

       output:
       set val(sample), file(sv), file(snv), file("*cluster_sorted_merged.bed") into beds
       

       shell:
       '''
       close_bp=!{params.close_value}
       closer_bp=!{params.closer_value}
       bp_per_kb=1000
       close=$((close_bp / bp_per_kb))
       closer=$((closer_bp / bp_per_kb))
    
       bcftools view -f 'PASS' !{sv} -Oz > !{sample}.sv.filt.vcf.gz
       Rscript !{baseDir}/vcf_to_bed.R --VCF !{sample}.sv.filt.vcf.gz --close !{params.close_value} --closer !{params.closer_value}
       sort -k1,1 -k2,2n !{sample}_0_${close}kb_cluster.bed > !{sample}_cluster_sorted.bed
       bedtools merge -i !{sample}_cluster_sorted.bed > !{sample}_cluster_sorted_merged.bed
             
       '''
  }


process make_vcfs {
    
    publishDir params.output_folder+"/SNV_clusters_VCFs/", mode: 'move', pattern: '*.snv.vcf*'    


    input:
    set val(sample), file(sv), file(snv),  file("*cluster_sorted_merged.bed") from beds
    file fasta_ref
    
    output:
    set val(sample), file(sv), file(snv), file("*snv*") into vcfs
    
    shell:
    '''
    close_bp=!{params.close_value}
    closer_bp=!{params.closer_value}
    bp_per_kb=1000
    close=$((close_bp / bp_per_kb))
    closer=$((closer_bp / bp_per_kb))
    echo $close $closer
   
    bcftools view -f 'PASS' !{snv} -Oz > !{sample}.filt.vcf.gz
    tabix -p vcf !{sample}.filt.vcf.gz
    
    bcftools view -f PASS --types snps --regions-file cluster_sorted_merged.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_cluster.snv.vcf.gz
    gunzip *.snv.vcf.gz
    '''
    
}   