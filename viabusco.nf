#!/usr/bin/env nextflow

if(!workflow.profile.contains('BUSCOs')) {
  println("This workflow requires -profile BUSCOs")
  exit 1
}

//INPUT PARAMS
trialLines = params.trialLines

import static groovy.json.JsonOutput.*
/*
  Generic method for extracting a string tag or a file basename from a metadata map
 */
 def getTagFromMeta(meta, delim = '_') {
  return (meta.species+delim+meta.version+(trialLines == null ? "" : delim+trialLines+delim+"trialLines"))
}

/*
  Generic method for extracting a string tag or a file basename from a metadata map
 */
def getAnnotationTagFromMeta(meta, delim = '_') {
  return meta.species+delim+meta.version+(meta.containsKey("annotation") ? delim+meta.annotation : "")+(trialLines == null ? "" : delim+trialLines+delim+"trialLines")
}


def helpMessage() {
  log.info"""
  Usage:

  nextflow run plantinformatics/pretzel-input-generator

  Default params:
  """.stripIndent()
  println(prettyPrint(toJson(params)))
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}


//ARRANGE INPUTS FOR PROCESSES
referencesLocal = Channel.create()
referencesRemote = Channel.create()
params.assemblies.each { assembly ->
  //sanitise any trailing spaces in strings
  // assembly.each {k,v ->  assembly[k] = v.trim()}
  //Abbreviate Genus_species name to G_species
  assembly.species = (assembly.species =~ /^./)[0]+(assembly.species =~ /_.*$/)[0]
  if(assembly.containsKey("fasta")) {
    if((assembly.fasta).matches("^(https?|ftp)://.*\$")) {
      referencesRemote << assembly
    } else {
      referencesLocal << [assembly,file(assembly.fasta)]
    }
  }
}
referencesRemote.close()
referencesLocal.close()


process fetchRemoteReference {
  // storeDir {executor == 'awsbatch' ? "${params.outdir}/downloaded" : "downloaded"}
  // scratch false

  tag{meta.subMap(['species','version'])}
  label 'download'


  input:
    val(meta) from referencesRemote

  output:
    set val(meta), file("${basename}.fasta") into referencesRemoteFasta

  script:
    basename=getTagFromMeta(meta)
    //DECOMPRESS?
    cmd = (meta.fasta).matches("^.*\\.gz\$") ?  "| gunzip --stdout " :  " "
    //TRIAL RUN? ONLY TAKE FIRST n LINES
    cmd += trialLines != null ? "| head -n ${trialLines}" : ""
    """
    curl -L ${meta.fasta} ${cmd} > ${basename}.fasta
    """
}

//Mix local and remote references then connect o multiple channels
// referencesRemoteFasta.mix(referencesLocal).into{ references }

process indexReference() {
  tag{meta.subMap(['species','version'])}
  label 'samtools'

  input:
    set val(meta), file(fasta) from referencesRemoteFasta.mix(referencesLocal)

  output:
    set val(meta), file(fasta), file("${fasta}.fai") into indexedReferences1, indexedReferences2

  """
  samtools faidx ${fasta}
  """
}


lineage = params.busco.lineage

process fetchAndPrepBuscoData {
  scratch false
  tag{ "${fname}" }

  input:
    val params.busco.lineage

  output:
    file (lineage.substring(lineage.lastIndexOf('/') + 1).replaceAll('.tar.gz','')) into lineageChannel

  script:
    fname = lineage.substring(lineage.lastIndexOf('/') + 1);
    //DOWNLOAD?
    fetch = lineage.matches("^(https?|ftp)://.*\$") ? "wget --no-check-certificate" : "ln -s"
    //DECOMPRESS?
    unzip = lineage.matches("^.*\\.tar\\.gz\$") ?  "tar xzvf ${fname}" :  " "
    """
    ${fetch} ${lineage}
    ${unzip}
    """
}

process runBUSCO {
  label 'BUSCO'
  label 'summary'

  tag{outmeta.subMap(['species','version','lineage'])}

  input:
    set val(meta), file(fasta), file("${fasta}.fai"), file(lineage) from indexedReferences2.combine(lineageChannel)

  output:
    set val(outmeta), file("run_${basename}/full_table_${basename}.tsv"), file("${fasta}.fai") into computedBUSCOs
    file ("run_${basename}/short_summary_${basename}.txt")

  script:
    outmeta = meta.clone()
    outmeta.lineage = lineage.name
    basename=getTagFromMeta(meta)
    //BUSCO wants to write to ${AUGUSTUS_CONFIG_PATH} which may be in a read-only container!
    """
    cp -r \${AUGUSTUS_CONFIG_PATH} augustus_config
    export AUGUSTUS_CONFIG_PATH=augustus_config
    run_BUSCO.py -i ${fasta} -o ${basename} --lineage_path ${lineage} --mode genome --cpu ${task.cpus} --species ${params.busco.augustusSpecies} --tarzip
    rm -r augustus_config
    """
}

process generateStandaloneJSONfromBUSCOs {
  tag{tag}
  label 'json'
  label 'groovy'
  errorStrategy 'ignore'

  input:
    set val(meta), file(tsv), file(fai) from computedBUSCOs

  output:
    file "*.json.gz" into featuresJSON
    // file "*.json.gz" into featuresJSON

  script:
    tag=getAnnotationTagFromMeta(meta)
    genome=getTagFromMeta(meta)
    shortName = (meta.containsKey("shortName") ? meta.shortName : "")
    shortName +=(meta.containsKey("annotation") ? "_"+meta.annotation : "") //only for cases where multiple annotations per genome
    template 'BUSCO_2_standaloneJSON.groovy'
}
