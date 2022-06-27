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
log.info "  AID mutations "
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
    log.info 'nextflow run main.nf '
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --input_file                   FILE           Input .csv file containing 2 columns: sample name and SNV vcf path.'
    log.info '    --output_folder                FOLDER         Output folder.'
    log.info '    --Off_targets_bed              FILE           Bed containing the coordinates to off targets of AID.'
    log.info '    --Off_flanking_bed              FILE          Bed containing the coordinates to off target flanking regions.'
    log.info '    --On_targets_bed               FILE           Bed containing the coordinates to correct/on targets of AID.'
    log.info '    --On_flaking_bed              FILE            Bed containing the coordinates to on target flanking regions.'
    log.info '    --fasta_ref                    FILE           Fasta reference file.'
    log.info ''
    log.info 'Flags:'
    log.info '    --help                                        Display this message'
    log.info ''
    exit 0
}

params.input_file = "/g/strcombio/fsupek_cancer1/SV_clusters_project/Hartwig_all_samples.csv"
params.Off_targets_bed = "/g/strcombio/fsupek_cancer1/SV_clusters_project/AID_beds/off_crg75_nochr.bed"
params.On_targets_bed = "/g/strcombio/fsupek_cancer1/SV_clusters_project/AID_beds/on_crg75_nochr.bed"
params.Off_flanking_bed = "/g/strcombio/fsupek_cancer1/SV_clusters_project/AID_beds/flanking_off_crg75_nochr_excl.bed"
params.On_flanking_bed = "/g/strcombio/fsupek_cancer1/SV_clusters_project/AID_beds/flanking_on_crg75_nochr_excl.bed"
params.output_folder = "/g/strcombio/fsupek_cancer1/SV_clusters_project/AID_mutations"
params.fasta_ref = "/g/strcombio/fsupek_cancer1/SV_clusters_project/hg19.fasta"

fasta_ref=file(params.fasta_ref)
Off_targets_bed=file(params.Off_targets_bed)
On_targets_bed=file(params.On_targets_bed)
Off_flanking_bed=file(params.Off_flanking_bed)
On_flanking_bed=file(params.On_flanking_bed)

pairs_list = Channel.fromPath(params.input_file, checkIfExists: true).splitCsv(header: true, sep: '\t', strip: true)
                   .map{ row -> [ row.sample, file(row.snv) ] }.view()
 

process get_vcfs {

       publishDir params.output_folder, mode: 'move', pattern: '*snv.vcf*'
       tag {sample}

       input:
       set val(sample), file(snv) from pairs_list
       file fasta_ref
       file Off_targets_bed
       file On_targets_bed
       file Off_flanking_bed
       file On_flanking_bed
    
       output:
       set val(sample), file(snv), file("*snv.vcf*") into vcfs

       shell:
       '''
       bcftools view -f 'PASS' !{snv} -Oz > !{sample}.filt.vcf.gz
       tabix -p vcf !{sample}.filt.vcf.gz
       
       bcftools view --types snps --regions-file !{Off_targets_bed} !{sample}.filt.vcf.gz | bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_Offtarget_region.snv.vcf.gz
       bcftools view --types snps --regions-file !{On_targets_bed} !{sample}.filt.vcf.gz | bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_Ontarget_region.snv.vcf.gz
       bcftools view --types snps --regions-file !{Off_flanking_bed} !{sample}.filt.vcf.gz | bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_Offtarget_flanking.snv.vcf.gz
       bcftools view --types snps --regions-file !{On_flanking_bed} !{sample}.filt.vcf.gz | bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_Ontarget_flanking.snv.vcf.gz
       gunzip *.snv.vcf.gz       
       '''
}
