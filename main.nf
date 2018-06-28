#!/usr/bin/env nextflow
echo true

//INPUT PARAMS
eprelease = params.eprelease

//LIST ELEMS ASSUMED TO MATCH
assemblySpeciesList = params.assemblySpecies.tokenize(",")
assemblyVersionList = params.assemblyVersion.tokenize(",")

//STATIC(?)
urlprefix = params.urlprefix
pepsuffix = params.pepsuffix
idxsuffix = params.idxsuffix

/*
* Download peptide seqs and assembly index files from Ensembl plants
*/
process fetchRemoteData {
  tag{species+" "+version}
  
  input:
    val species from assemblySpeciesList
    val version from assemblyVersionList

  output:
    set val(species), val(version), file(idx) into remoteIndices
    set val(species), val(version), file(pep) into remotePepSeqs4Features
    
 
  script:

    idxurl=urlprefix+eprelease+"/fasta/"+species.toLowerCase()+"/dna_index/"+species+"."+version+idxsuffix
    pepurl=urlprefix+eprelease+"/fasta/"+species.toLowerCase()+"/pep/"+species+"."+version+pepsuffix
    """
    curl $idxurl > idx
    curl $pepurl | gunzip --stdout  > pep
    """
}

/*
* Generate genome blocks definitions JSON for pretzel 
*/
process generateGenomeBlocksJSON {
  tag{tag}
  label 'json'
  //publishDir "${params.outdir}/JSON"

  input:
    set val(species), val(version), file(idx) from remoteIndices

  output:
     file "*.json"

  script: 
    tag=species+"_"+version
    """
    awk '\$1 ~/^(chr|[0-9])/' "${idx}" \
    | faidx2json.awk -vname="${tag}" > "${tag}"_genome.json
    """
}

/*
* Generate for pretzel JSON aliases linking features between chromosomes/genomes
*/
process generateFeaturesJSON {
  tag{tag}
  label 'json'

  input:
    set val(species), val(version), file(pep) from remotePepSeqs4Features
  
   output:
    file "*.json"

  script:    
    tag=species+"_"+version
    """
    awk '\$1 ~ /^>/' "${pep}"  | sort -k3,3V | pep2featuresJSON.awk -vname="${tag}" > "${tag}"_annotation.json
    """
}

// /*
// * Generate for pretzel JSON aliases linking features between chromosomes/genomes
// */
// process generateAliasesJSON {
//   input:
//     set val(species), val(version), file(pep) from remotePepSeqs

// }

// process fetchReads {

//   input: 
//     val reads1url
//     val reads2url

//   output:
//     set val(longtag), val(nametag),file("r1.gz"), file("r2.gz") into FASTQ, hisat2FASTQ, kangaFASTQ

//   script:
//     nametag = "tmpTAG"
//     longtag = ["name":"real", "nreads":"10000", "seqerr":"NaN", "rep":"na", "format":"fq"]
//     """
//     curl ${reads1url} | gunzip --stdout | head -n 40000 | pigz --fast > r1.gz
//     curl ${reads2url} | gunzip --stdout | head -n 40000 | pigz --fast > r2.gz
//     """

// }

