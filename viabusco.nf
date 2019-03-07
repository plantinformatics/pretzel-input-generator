#!/usr/bin/env nextflow

if(!workflow.profile.contains('BUSCO')) {
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

// /*
// * Generate genome blocks definitions JSON for pretzel
// */
// process generateGenomeBlocksJSON {
//   tag{tag}
//   label 'json'
//   label 'groovy'

//   input:
//     set val(meta), file(fasta), file(idx) from indexedReferences1

//   output:
//     file "*.json" into genomeBlocksJSON

//   script:
//     tag=getTagFromMeta(meta)
//     """
//     #!/usr/bin/env groovy

//     import static groovy.json.JsonOutput.*
//     idx = new File('${idx}').text
//     out = new File('${tag}_genome.json')
//     def genome = [:]
//     genome.name = "${tag}"
//     genome.meta = [:]
//     if("${meta.shortName}" != "null") {
//       genome.meta << ["shortName" : "${meta.shortName}"]
//     }
//     if("${meta.source}" != "null") {
//       genome.meta << ["source" : "${meta.source}"]
//     }
//     if("${meta.version}" != "null") {
//       genome.meta << ["version" : "${meta.version}"]
//     }
//     if("${meta.citation}" != "null") {
//       genome.meta << ["citation" : "${meta.citation}"]
//     }
//     genome.blocks = []
//     idx.eachLine { line ->
//       if(line.toLowerCase() =~ /^(chr|[0-9]|x|y|hic_scaffold)/ ) {
//         toks = line.split('\t')
//         genome.blocks += [ "scope": toks[0].replaceFirst("^(C|c)(H|h)(R|r)[_]?",""), "featureType": "linear", "range": [1, toks[1].toInteger()] ]
//       }
//     }
//     out.text = prettyPrint(toJson(genome))
//     """
// }


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
  label 'tsv'

  tag{outmeta.subMap(['species','version','lineage'])}

  input:
    set val(meta), file(fasta), file("${fasta}.fai"), file(lineage) from indexedReferences2.combine(lineageChannel)

  output:
    set val(outmeta), file("run_${basename}/full_table_${basename}.tsv"), file("${fasta}.fai") into computedBUSCOs

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


// process generateFeaturesJSONfromBUSCOs {
//   tag{tag}
//   label 'json'
//   label 'groovy'

//   input:
//     set val(meta), file(tsv), file(fai) from computedBUSCOs

//   output:
//     file "*.json.gz" into featuresJSON

//   script:
//     tag=getAnnotationTagFromMeta(meta)
//     genome=getTagFromMeta(meta)
//     shortName = (meta.containsKey("shortName") ? meta.shortName+"_genes" : "")
//     shortName +=(meta.containsKey("annotation") ? "_"+meta.annotation : "") //only for cases where multiple annotations per genome
//     """
//     #!/usr/bin/env groovy

//     import static groovy.json.JsonOutput.*
//     buscos = new File('${tsv}').text
//     idx = new File('${fai}').text
//     out = new File('${tag}_annotation.json')
//     //RECORD chromosome/scaffold sizes
//     lengths. = [:]
//     idx.eachLine { line ->
//       if(line.toLowerCase() =~ /^(chr|[0-9]|x|y|hic_scaffold)/ ) {
//         toks = line.split('\t')
//         lengths.(toks[0].replaceFirst("^(C|c)(H|h)(R|r)[_]?","") = toks[1].toInteger()
//       }
//     }
//     def annotation = [:]
//     annotation.meta = [:]
//     if("${meta.shortName}" != "null") {
//       annotation.meta << ["shortName" : "${shortName}"]
//     }
//     if("${meta.source}" != "null") {
//       annotation.meta << ["source" : "${meta.source}"]
//     }
//     if("${meta.version}" != "null") {
//       annotation.meta << ["version" : "${meta.version}"]
//     }
//     if("${meta.citation}" != "null") {
//       annotation.meta << ["citation" : "${meta.citation}"]
//     }
//     annotation.name = "${tag}_genes"
//     annotation.namespace = "${genome}:${tag}_annotation"
//     //annotation.namespace = "${meta.lineage}"
//     annotation.parent = "${genome}"
//     annotation.blocks = []
//     TreeMap scope = [:] //keep keys sorted as the corresponding blocks get displayed in order in pretzel
//     buscos.eachLine { line ->
//       if(!line.startsWith('#')) {
//         toks = line.split('\\t')
//         if(toks[1].matches("${params.busco.allowedStatus}") ) {
//           location = toks[2].split(":")
//           gene = toks[0]
//           key = toks[2]
//           if(!scope.containsKey(key)) {
//             scope << [(key) : []]
//           }
//           scope[key] << ["name" : gene, "range" : [ toks[3].toInteger(), toks[4].toInteger() ]]
//         }
//       }
//     }
//     //GROUP TOGETHER FEATURES FROM/IN SAME BLOCK
//     scope.each { k, features ->
//       if(features.size() >= ${params.busco.minPlaced}) {
//         current = [ "scope": k, "featureType": "linear", "range": [1, lengths[k]], "features": []]
//         features.each { feature ->
//           current.features << feature
//         }
//         annotation.blocks << current
//       }
//     }
//     out.text = prettyPrint(toJson(annotation))
//     'gzip ${tag}_annotation.json'.execute()
//     """
// }

process generateStandaloneJSONfromBUSCOs {
  tag{tag}
  label 'json'
  label 'groovy'
  errorStrategy 'ignore'

  input:
    set val(meta), file(tsv), file(fai) from computedBUSCOs

  output:
    file "*.json" into featuresJSON
    // file "*.json.gz" into featuresJSON

  script:
    tag=getAnnotationTagFromMeta(meta)
    genome=getTagFromMeta(meta)
    shortName = (meta.containsKey("shortName") ? meta.shortName : "")
    shortName +=(meta.containsKey("annotation") ? "_"+meta.annotation : "") //only for cases where multiple annotations per genome
    """
    #!/usr/bin/env groovy

    import static groovy.json.JsonOutput.*
    buscos = new File('${tsv}').text
    out = new File('${tag}_annotation.json')
    def annotation = [:]
    annotation.meta = [:]
    if("${meta.shortName}" != "null") {
      annotation.meta << ["shortName" : "${shortName}"]
    }
    if("${meta.source}" != "null") {
      annotation.meta << ["source" : "${meta.source}"]
    }
    if("${meta.version}" != "null") {
      annotation.meta << ["version" : "${meta.version}"]
    }
    if("${meta.citation}" != "null") {
      annotation.meta << ["citation" : "${meta.citation}"]
    }
    annotation.name = "${tag}"
    annotation.namespace = "${meta.lineage}"
    annotation.blocks = []
    TreeMap scope = [:] //keep keys sorted as the corresponding blocks get displayed in order in pretzel
    buscos.eachLine { line ->
      if(!line.startsWith('#')) {
        toks = line.split('\\t')
        if(toks[1].matches("${params.busco.allowedStatus}") ) {
          location = toks[2].split(":")
          gene = toks[0]
          key = toks[2]
          if(!scope.containsKey(key)) {
            scope << [(key) : []]
          }
          scope[key] << ["name" : gene, "range" : [ toks[3].toInteger(), toks[4].toInteger() ]]
        }
      }
    }
    //GROUP TOGETHER FEATURES FROM/IN SAME BLOCK
    scope.each { k, features ->
      if(features.size() >= ${params.busco.minPlaced}) {
        current = [ "scope": k, "featureType": "linear", "features": []]
        features.each { feature ->
          current.features << feature
        }
        annotation.blocks << current
      }
    }
    out.text = prettyPrint(toJson(annotation))
    // 'gzip ${tag}_annotation.json'.execute()
    """
}