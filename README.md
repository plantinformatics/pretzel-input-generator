# pretzel-input-generator

This is a [nextflow](https://www.nextflow.io) pipeline for generating input for [pretzel](https://github.com/plantinformatics/pretzel) from annotated and (mostly) contiguous genome assemblies. 

## Input



### Remote

The pipeline pulls data from [Ensembl plants](https://plants.ensembl.org/index.html), included species and assembly versions are specified in [conf/input.config](conf/input.config). 

* genome assembly index file 
* matching protein sequences 

### Local

The pipeline requires 

* a genome assembly index file
* gene annotations (currently GTF but will also accept GFF3)
* matching protein sequences

## Dependencies

* [nextflow](https://www.nextflow.io) 
* Either of:
  * Singularity
  * Docker
  * Required software installed (e.g. as module, in which case specify its name in [`conf/modules.config`](conf/modules.config)) 

When using Singularity or Docker, the required containers are specified in [`conf/containers.conf`](conf/containers.config)
 
<!-- [MMSeqs2](https://github.com/soedinglab/mmseqs2) -->


## Execution

We provide several execution profiles, "locally" may mean a designated server or an interactive session on a cluster. By appending  e.g. `-revision v0.1` to your command you can specify a release tag to run a revision that may be more likely to be working than the latest.

Run locally with docker

```
nextflow run rsuchecki/pretzel-input-generator -resume -profile docker
```

Run locally with singularity

```
nextflow run rsuchecki/pretzel-input-generator -resume -profile singularity
```

Dispatch on a SLURM cluster with singularity

```
nextflow run rsuchecki/pretzel-input-generator -resume -profile slurm,singularity,singularitymodule
```

Dispatch on a SLURM cluster with modules

```
nextflow run rsuchecki/pretzel-input-generator -resume -profile slurm,modules
```