#!/usr/bin/env nextflow

//INPUT PARAMS
trialLines = params.trialLines

import static groovy.json.JsonOutput.*
/*
  Generic method for extracting a string tag or a file basename from a metadata map
 */
 def getTagFromMeta(meta, delim = '_') {
  return (meta.species+delim+meta.version+(trialLines == null ? "" : delim+trialLines+delim+"trialLines"))
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
    set val(meta), file(fasta), file("${fasta}.fai") into indexedReferences

  """
  samtools faidx ${fasta}
  """
}

lineage = params.lineageBUSCO

process fetchAndPrepBuscoData {
  scratch false
  tag{ "${fname}" }

  input:
    val lineage

  output:
    file (lineage.substring(lineage.lastIndexOf('/') + 1).replaceAll('.tar.gz','')) into lineageChannel

  script:
    fname = lineage.substring(lineage.lastIndexOf('/') + 1);
    //DOWNLOAD?
    fetch = lineage.matches("^(https?|ftp)://.*\$") ? "wget" : "ln -s"
    //DECOMPRESS?
    unzip = lineage.matches("^.*\\.tar\\.gz\$") ?  "tar xzvf ${fname}" :  " "
    """
    ${fetch} ${lineage}
    ${unzip}
    """
}

process runBUSCO {
  label 'BUSCO'
  tag{outmeta.subMap(['species','version','lineage'])}

  input:
    set val(meta), file(fasta), file("${fasta}.fai"), file(lineage) from indexedReferences.combine(lineageChannel).first()

  // output:
  //   set val(outmeta),

  script:
    outmeta = meta.clone()
    outmeta.lineage = lineage
    basename=getTagFromMeta(meta)
    //BUSCO wants to write to ${AUGUSTUS_CONFIG_PATH} which may be in a read-only container!
    """
    cp -r \${AUGUSTUS_CONFIG_PATH} augustus_config
    export AUGUSTUS_CONFIG_PATH=augustus_config
    run_BUSCO.py -i ${fasta} -o ${basename} --lineage_path ${lineage} --mode genome --cpu ${task.cpus} --species ${params.augustusSpecies} --tarzip
    rm -r augustus_config
    """
}


