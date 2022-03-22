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
    log.info 'nextflow run main.nf --input_file list_SV_SNV.txt --distance 1'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --input_file                   STRING         Input .csv file containing 3 columns: sample name, SV vcf path and SNV vcf path.'
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

params.input_file = null
params.closer_value = null
params.close_value = null
params.output_folder = null 
params.fasta_ref = null
params.hg19 = NULL 


pairs_list = Channel.fromPath(params.input_file, checkIfExists: true).splitCsv(header: true, sep: '\t', strip: true)
                   .map{ row -> [ row.sample, file(row.sv), file(row.snv) ] }.view()
hg19 = file(params.hg19)


process make_beds {

       publishDir params.output_folder+"/SV_BEDs/", mode: 'copy', pattern: '*.bed'

       tag {sample}

       input:
       set val(sample), file(sv), file(snv) from pairs_list
       file hg19

       output:
       set val(sample), file(sv), file(snv) into beds
       set file("*unclustered.bed") file("*

       shell:
       '''
       bcftools view -f 'PASS' !{sv} -Oz > !{sample}_filt.vcf.gz
       Rscript  !{baseDir}/vcf_to_bed.R --VCF !{sample}_filt.vcf.gz --close !{params.close_value} --closer !{params.closer_value}
       bedtools complement -i !{sample}_0_!{params.close_value}.bed -g !{hg19} > !{sample}_unclustered.bed
       '''
  }


process make_vcfs {
    
    publishDir params.output_folder+"/SNV_clusters_VCFs/", mode: 'move', pattern: '*_clustered_*.vcf.gz'
    publishDir params.output_folder+"/SNV_clusters_VCFs/", mode: 'move', pattern: '*_unclustered_*.vcf.gz'

    
    tag {sample}


    
    
    
}   
    
    
