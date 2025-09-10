process hello {
  conda '/home/dmmalone/miniconda3/envs/NIWWt'

  script:
  """
  trim_galore -v
  """
}

workflow {
    hello()
}