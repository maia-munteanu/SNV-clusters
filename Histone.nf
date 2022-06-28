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
log.info "--------------------------------------------------------------------------"
log.info "  Histone mutations - extracting mutations from regions marked by histones"
log.info "--------------------------------------------------------------------------"
log.info "Copyright (C) IRB Barcelona"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------------------------"
log.info ""

params.help = null
if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info '  USAGE              '
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run AID.nf '
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --input_file                   FILE           Input .csv file containing 2 columns: sample name and SNV vcf path.'
    log.info '    --output_folder                FOLDER         Output folder.'
    log.info '    --K4_K27                       FILE           Intersected beds.'
    log.info '    --K4_K36                       FILE           Intersected beds.'
    log.info '    --K27_K36                      FILE           Intersected beds.'
    log.info '    --K4_K27_K36                   FILE           Intersected beds.'
    log.info '    --fasta_ref                    FILE           Fasta reference file.'
    log.info '    --CRG75                        FILE           Low mappability bed, here CRG75 where all regions with >1 occurance are excluded.'
    log.info ''
    log.info 'Flags:'
    log.info '    --help                                        Display this message'
    log.info ''
    exit 0
}

params.input_file = "/g/strcombio/fsupek_cancer1/SV_clusters_project/Hartwig_all_samples.csv"
params.output_folder = "/g/strcombio/fsupek_cancer1/SV_clusters_project/Histone_mutations"
params.fasta_ref = "/g/strcombio/fsupek_cancer1/SV_clusters_project/hg19.fasta"
params.CRG75 = "/home/mmunteanu/reference/CRG75_nochr.bed"
params.K4_K27 = "/g/strcombio/fsupek_cancer1/SV_clusters_project/Histone_beds/pooled8Solid_H3K4me3_H3K27ac.3of3.nochr.bed"
params.K4_K36 = "/g/strcombio/fsupek_cancer1/SV_clusters_project/Histone_beds/pooled8Solid_H3K4me3_H3K36me3.3of3.nochr.bed"
params.K27_K36 = "/g/strcombio/fsupek_cancer1/SV_clusters_project/Histone_beds/pooled8Solid_H3K27ac_H3K36me3.3of3.nochr.bed"
params.K4_K27_K36 = "/g/strcombio/fsupek_cancer1/SV_clusters_project/Histone_beds/pooled8Solid_H3K4me3_H3K27ac_H3K36me3.3of3.nochr.bed"

fasta_ref=file(params.fasta_ref)
CRG75=file(params.CRG75)
K4_K27=file(params.K4_K27)
K4_K36=file(params.K4_K36)
K27_K36=file(params.K27_K36)
K4_K27_K36=file(params.K4_K27_K36)


pairs_list = Channel.fromPath(params.input_file, checkIfExists: true).splitCsv(header: true, sep: '\t', strip: true)
                   .map{ row -> [ row.sample, file(row.snv) ] }.view()
 

process get_vcfs {

       publishDir params.output_folder, mode: 'move', pattern: '*snv.vcf*'
       tag {sample}

       input:
       set val(sample), file(snv) from pairs_list
       file fasta_ref
       file CRG75
       file K4_K27
       file K4_K36
       file K27_K36
       file K4_K27_K36
    
       output:
       set val(sample), file(snv), file("*snv.vcf*") into vcfs

       shell:
       '''
       tabix -p vcf !{snv}
       bcftools view -f 'PASS' --regions-file !{CRG75} !{snv} | bcftools sort -Oz > !{sample}.filt.vcf.gz
       tabix -p vcf !{sample}.filt.vcf.gz
       
       bcftools view --types snps --regions-file !{K4_K27} !{sample}.filt.vcf.gz | bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_H3K4me3_H3K27ac.3of3.snv.vcf.gz
       bcftools view --types snps --regions-file !{K4_K36} !{sample}.filt.vcf.gz | bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_H3K4me3_H3K36me3.3of3.snv.vcf.gz
       bcftools view --types snps --regions-file !{K27_K36} !{sample}.filt.vcf.gz | bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_H3K27ac_H3K36me3.3of3.snv.vcf.gz
       bcftools view --types snps --regions-file !{K4_K27_K36} !{sample}.filt.vcf.gz | bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_H3K4me3_H3K27ac_H3K36me3.3of3.snv.vcf.gz
       gunzip *.snv.vcf.gz       
       '''
}
