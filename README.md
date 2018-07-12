# pretzel-input-generator

This is a [nextflow](https://www.nextflow.io) pipeline for generating input for [pretzel](https://github.com/plantinformatics/pretzel) from annotated and (mostly) contiguous genome assemblies. 

## Input

Input files are specified in [conf/input.config](conf/input.config). This can be supplemented/replaced by JSON/YAML formatted input spec.  

### Remote

The pipeline pulls data from [Ensembl plants](https://plants.ensembl.org/index.html), included species and assembly versions are specified in [conf/input.config](conf/input.config). 
For each of the datsets the pipeline downloads:

* genome assembly index file 
* matching protein sequences 

### Local

The pipeline requires 

* a genome assembly index file
* gene annotations (either GTF or GFF3)
* matching protein sequences (presumably for representative isoform) 

Currently, some datasets are specified and to run this pipeline without the specified files, use `--localAssembly NA` at execution or modify your local copy of [conf/input.config](conf/input.config)

## Dependencies

* [nextflow](https://www.nextflow.io) 
* Either of:
  * Singularity
  * Docker
  * Required software installed (e.g. as module, in which case specify its name in [`conf/modules.config`](conf/modules.config)) 

When using Singularity or Docker, the required containers are specified in [`conf/containers.conf`](conf/containers.config)
 
<!-- [MMSeqs2](https://github.com/soedinglab/mmseqs2) -->


## Execution

We provide several execution profiles, "locally" may mean a designated server or an interactive session on a cluster. By appending  e.g. `-revision v0.2` to your command you can specify a release tag to run a revision that may be more likely to be working than the latest.

Run locally with docker

```
nextflow run rsuchecki/pretzel-input-generator -resume -profile docker --localAssembly NA
```

Run locally with singularity

```
nextflow run rsuchecki/pretzel-input-generator -resume -profile singularity --localAssembly NA
```

Dispatch on a SLURM cluster with singularity

```
nextflow run rsuchecki/pretzel-input-generator -resume -profile slurm,singularity,singularitymodule --localAssembly NA
```

Dispatch on a SLURM cluster with modules

```
nextflow run rsuchecki/pretzel-input-generator -resume -profile slurm,modules --localAssembly NA
```

