params.refFa = '~/ansel_project1/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa'
params.refGtf = '~/ansel_project1/reference/Homo_sapiens.GRCh38.91.gtf'
params.reads = '~/ansel_project1/fastq/SGNex_H9_directRNA_replicate3_run1.fastq.gz'
params.outdir = '~/ansel_project1/results'

process MINIMAP2_ALIGN {
  input:
    path refFa
    path reads
  output:
    path "aligned_reads.sam"

  """
    minimap2 -ax splice -uf -k14 $refFa $reads > aligned_reads.sam
  """
}

process SAM_TO_BAM {
  publishDir params.outdir, mode: 'copy'

  input:
    path reads_sam
  output:
    path "aligned_reads.bam"

  """
    samtools view -b $reads_sam > aligned_reads.bam
  """
}

process BAMBU {
  publishDir params.outdir, mode: 'copy'

  input:
    path refFa
    path refGtf
    path reads_bam
  output:
    path "counts_transcript.txt"
    path "counts_gene.txt"
    path "extended_annotations.gtf"

    """
    #!/usr/bin/env Rscript --vanilla
    library(bambu)
    annotations <- prepareAnnotations("$refGtf")
    se     <- bambu(reads = "$reads_bam",
                    annotations = annotations,
                    genome = "$refFa",
                    NDR=1,
                    ncore = 1)
    writeBambuOutput(se, path = "./")
    """
}

workflow {
  MINIMAP2_ALIGN(params.refFa, params.reads)
  SAM_TO_BAM(MINIMAP2_ALIGN.out)
  BAMBU(params.refFa, params.refGtf, SAM_TO_BAM.out)
}
