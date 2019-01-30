[![Latest GitHub tag](https://img.shields.io/github/tag/plantinformatics/pretzel-input-generator.svg?label=latest%20release&logo=github&style=for-the-badge)](https://github.com/plantinformatics/pretzel-input-generator/releases)

# Pipeline overview

`pretzel-input-generator` is a [nextflow](https://www.nextflow.io) pipeline for generating input for [pretzel](https://github.com/plantinformatics/pretzel) from annotated and (mostly) contiguous genome assemblies. The pipeline requires approximately 1 cpu-day, but as many processes can run independently, the real run-time is much shorter if suitable compute resources are available.


<!-- TOC -->

- [Pipeline overview](#pipeline-overview)
- [Default pipeline](#default-pipeline)
  - [Quick start](#quick-start)
  - [Input](#input)
    - [Data sources](#data-sources)
      - [Remote](#remote)
      - [Local](#local)
  - [Dependencies](#dependencies)
  - [Execution](#execution)
  - [Output](#output)
- [BUSCO-based pipeline](#busco-based-pipeline)
  - [Quick-ish start](#quick-ish-start)
  - [Output](#output-1)

<!-- /TOC -->

# Default pipeline

Designed for EnsemblPlants and similarily formatted data.

![doc/dag.png](doc/dag.png)


## Quick start

Requires [nextflow](https://www.nextflow.io) and [Singularity](http://singularity.lbl.gov)

```
nextflow run plantinformatics/pretzel-input-generator \
-revision v1.1 -profile EP,singularity --localAssembly NA
```

This will pull and process data sets from [Ensembl plants](https://plants.ensembl.org) specified in [`conf/triticeae.config`](conf/input.triticeae#L9-L29)

## Input

Input files are specified in [conf/triticeae.config](conf/triticeae.config). This can be supplemented/replaced by JSON/YAML formatted input spec.

### Data sources

Currently all input data comes from the following sources:

* [Ensembl plants](https://plants.ensembl.org) - multiple datasets as specified in [`conf/input.config`](conf/input.config#L9-L29)
* [International Wheat Genome Sequencing Consortium](https://www.wheatgenome.org/)
  * [Triticum aestivum (Chinese Spring) IWGSC RefSeq v1.0 assembly](https://wheat-urgi.versailles.inra.fr/Seq-Repository/Assemblies)
* [The wild emmer wheat sequencing consortium (WEWseq)](http://wewseq.wixsite.com/consortium)
  * Zavitan assembly downloaded from [GrainGenes](https://wheat.pw.usda.gov/GG3/wildemmer)
* [European Nucleotide Archive](https://www.ebi.ac.uk/ena)
  * [Assembly of chromosome 2D of *Triticum aestivum* line CH Campala *Lr22a*](https://www.ebi.ac.uk/ena/data/view/LS480641)
  * [Assembly of *Triticum urartu* ](https://www.ebi.ac.uk/ena/data/view/GCA_003073215)
    * Annotation downloaded from [MBKBase](http://www.mbkbase.org/Tu/)
  * [Assembly of *Aegilops tauschii* ](https://www.ebi.ac.uk/ena/data/view/GCA_002575655.1)
    * Annotation downloaded from [http://aegilops.wheat.ucdavis.edu/ATGSP/annotation/](http://aegilops.wheat.ucdavis.edu/ATGSP/annotation/)

#### Remote

The pipeline pulls data from [Ensembl plants](https://plants.ensembl.org), included species and assembly versions are specified in [conf/triticeae.config](conf/triticeae.config).
For each of the datsets the pipeline downloads:

* genome assembly index file
* protein sequences

#### Local

The pipeline requires

* a genome assembly index file - all we need are lengths of pseudo-chromosomes so a two-column `.tsv` file with chromosome names and their lengths will suffice
* gene annotations (either GTF or GFF3)
* matching protein sequences (presumably for representative isoform)

If GTF/GFF3 is not available, the protein sequences FASTA id and description lines must be formatted to contain information as per the following example:

```
>AT1G24405.1 pep chromosome:TAIR10:1:8654945:8655662:1 gene:AT1G24405
```

This follows how protein sequences are annotated on Ensembl plants, but we do not **currently** use all the information in the description line, the complete version of which is:

```
>AT1G24405.1 pep chromosome:TAIR10:1:8654945:8655662:1 gene:AT1G24405 transcript:AT1G24405.1 gene_biotype:protein_coding transcript_biotype:protein_coding description:F21J9.7 [Source:UniProtKB/TrEMBL;Acc:Q9FYM2]
```

To run this pipeline without the local input files, use `--localAssembly NA` at execution or modify your local copy of [conf/triticeae.config](conf/triticeae.config)

Wherever possible the local assembly files are used as input for the pipeline in their original form - as downloaded from their respective sources. This is however not always possible due to inconsistencies in formatting and varying levels of adherence to standards and conventions. We try to capture additional steps needed to prepare these input data sets for the inclusion in this pipeline in [doc/format_local.md](doc/format_local.md).

## Dependencies

* [nextflow](https://www.nextflow.io)
* **Either** of the following:
  * [Singularity](http://singularity.lbl.gov)
  * [Docker](http://singularity.lbl.gov)
  * Required software installed. In addition to standard linux tools, these include:
    * [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
    * [MMSeqs2](https://github.com/soedinglab/mmseqs2)

When using Singularity or Docker, the required containers are specified in [`conf/containers.conf`](conf/containers.config)


## Execution

We provide several execution profiles, "locally" may mean a designated server or an interactive session on a cluster. By appending  e.g. `-revision v1.1` to your command you can specify a release tag to run a specific revision. When re-running the pipeline after errors or changes use `-resume` to ensure only the necessary processes are re-run.

Run locally with docker

```
nextflow run plantinformatics/pretzel-input-generator \
-profile EP,docker --localAssembly NA
```

Run locally with singularity

```
nextflow run plantinformatics/pretzel-input-generator \
-profile EP,singularity --localAssembly NA
```

Dispatch on a SLURM cluster with singularity

```
nextflow run plantinformatics/pretzel-input-generator \
-profile EP,slurm,singularity,singularitymodule --localAssembly NA
```

## Output

All generated JSON files generated by the pipeline are output to `results/JSON`.

* For each of the input genome assemblies, these include:
  * `*_genome.json` - dataset (genome) definitions specifying outer coordinates of blocks (chromosomes)
  * `*_annotation.json.gz` - specifications of coordinates of features (genes) within blocks
* In addition, for each (lexicographically ordered) pair of genome assemblies, the pipeline generates:
  * `*_aliases.json.gz` which specify links between features between the two genomes.

The output files (hopefully) conform to the requirements of [pretzel data structure](https://github.com/plantinformatics/pretzel-data).


The `results/flowinfo` directory contains summaries of pipeline execution and `results/downloads` includes the files downloaded from Ensembl plants.

```
results
├── downloads
├── flowinfo
└── JSON
```

To upload the generated data to your instance of pretzel, follow [these instructions](doc/upload.md).


# BUSCO-based pipeline

This approach is much simpler and at the same time computationally intensive.
Its main advantage is that it dos not require gene annotations, all that is required is a set of genome assemblies.

![doc/dag-busco.png](doc/dag-busco.png)

## Quick-ish start

```
nextflow run plantinformatics/pretzel-input-generator/viabusco.nf \
-profile BUSCOs,singularity
```

This will pull and process data sets from [DNA Zoo](https://www.dnazoo.org/) specified in [`conf/dna_zoo_felidae.config`](conf/dna_zoo_felidae.config), consuming around 280 CPU hours and given sufficient resources should complete in a day or so.

## Output

In comparison with the main pipeline the output lacks `*_aliases.json.gz` as features on different genomes are implicitly connected by BUSCOs identifiers.