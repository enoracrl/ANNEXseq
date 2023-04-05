process BAMBU_SPLIT_RESULTS {
  input:
  file extended_annotation

  output:
  path 'novel_genes.gtf', emit: novel_genes
  path 'novel_isoforms.gtf', emit: novel_isoforms

  shell:
  '''
  grep "BambuTx" !{extended_annotation} | awk '$3=="exon"' > novel.gtf
  grep "BambuGene" novel.gtf > novel_genes.gtf
  grep -v "BambuGene" novel.gtf > novel_isoforms.gtf
  '''
}
