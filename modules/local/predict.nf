process PREDICT {
  time '24.h'
  cpus 6
  memory '84 GB'

  conda (params.enable_conda ? "$baseDir/environment.yml" : null)
  container "ghcr.io/igdrion/annexa:${workflow.revision? workflow.revision: "dev"}"
  publishDir "$params.outdir/transforkmers", mode: 'copy'

  input:
  file tss_sequences
  path tokenizer
  path model

  output:
  path "output.csv", emit: tss_prob

  """
  transforkmers predict \
    --model_path_or_name ${model} \
    --tokenizer ${tokenizer} \
    --input ${tss_sequences} \
    --quantize-model \
    --output . \
    --per_device_eval_batch_size 64
  """
}