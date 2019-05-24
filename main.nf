#!/usr/bin/env nextflow

if(!workflow.profile.contains('EP')) {
  println("This workflow requires -profile EP to be specified")
  exit 1
}

//INPUT PARAMS
trialLines = params.trialLines
// eprelease = params.eprelease

import static groovy.json.JsonOutput.*

def helpMessage() {
  log.info"""
  =============================================================
  plantinformatics/pretzel-input-generator  ~  version ${params.version}
  =============================================================
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

//STATIC(?) ENSEMBL URLs
urlprefix = params.urlprefix
pepsuffix = params.pepsuffix
idxsuffix = params.idxsuffix

// def noKeyOrIsMap

//LOCAL INPUTS
localInput = Channel.create()
localIndices = Channel.create()
if(params.localAssembly != "NA") {
  params.localAssembly.each {
    // println(prettyPrint(toJson(it)))
    // for (key in ["gtf","gff3"]) {
    key = 'gtfgff3'
    // if(it.containsKey(key)) {
    //IF MORE THAN ONE ANNOTATION PER GENOME
    if(it.pep instanceof Map) { // && it.containsKey(key) && it.get(key) instanceof Map) {
      it.pep.each {id, pep ->
      // println(id+" -> "+pep)
        clone = it.clone()
        clone.annotation = id
        gtfgff3 = it.containsKey(key) ? file(it.get(key).get(id)) : null
          localInput << [clone,gtfgff3,file(pep)]
      }
    } else { //if (!(it.pep instanceof Map) && !(it.get(key) instanceof Map)){
      gtfgff3 = it.containsKey(key) ? file(it.get(key)) : null
        localInput << [it,gtfgff3,file(it.pep)]
      // } else {
      //   exitOnInputMismatch(it)
    }
    //ALL SHOULD HAVE AN INDEX
    if(it.containsKey("idx")) {
      localIndices << [it,file(it.idx)]
    }
  }
}
localInput.close()
localIndices.close()

// def exitOnInputMismatch(data) {
//   println("Malformed input. Expecting number of pep and gtf/gff3 inputs to match for a data set.")
//   println("Offending data set: ")
//   println(data)
//   println("Terminating.")
//   System.exit(1)
// }


/*
  Generic method for extracting a string tag or a file basename from a metadata map
 */
def getAnnotationTagFromMeta(meta, delim = '_') {
  return meta.species+delim+meta.version+(meta.containsKey("annotation") ? delim+meta.annotation : "")+(trialLines == null ? "" : delim+trialLines+delim+"trialLines")
}

/*
  Generic method for extracting a string tag or a file basename from a metadata map
 */
def getDatasetTagFromMeta(meta, delim = '_') {
  return meta.species+delim+meta.version+(trialLines == null ? "" : delim+trialLines+delim+"trialLines")
}

/*
* Download peptide seqs and assembly index files from Ensembl plants
*/
process fetchRemoteDataFromEnsemblPlants {
  tag{meta}
  label 'download'

  input:
    set val(species), val(version), val(shortName), val(eprelease) from Channel.from(params.remoteAssembly)

  output:
    set val(meta), file("${basename}.idx") into remoteIndices
    set val(meta), file("${basename}.pep") into remotePepSeqs

  script:
    meta=["species":species, "version":version, "source": "https://plants.ensembl.org/"+species, "release": eprelease, "shortName": shortName]
    basename=getDatasetTagFromMeta(meta)
    idxurl=urlprefix+eprelease+"/fasta/"+species.toLowerCase()+"/dna_index/"+species+"."+version+idxsuffix
    pepurl=urlprefix+eprelease+"/fasta/"+species.toLowerCase()+"/pep/"+species+"."+version+pepsuffix
    if(trialLines == null) {
      """
      curl $idxurl > ${basename}.idx
      curl $pepurl | gunzip --stdout > ${basename}.pep
      """
    } else {
      """
      curl $idxurl > ${basename}.idx
      curl $pepurl | gunzip --stdout | head -n ${trialLines} > ${basename}.pep
      """
    }
}

/*
* Generate genome blocks definitions JSON for pretzel
*/
process generateGenomeBlocksJSON {
  tag{tag}
  label 'json'
  label 'groovy'

  input:
    set val(meta), file(idx) from localIndices.mix(remoteIndices)

  output:
    file "*.json" into genomeBlocksJSON

  script:
    tag=getDatasetTagFromMeta(meta)
    """
    #!/usr/bin/env groovy

    import static groovy.json.JsonOutput.*
    idx = new File('${idx}').text
    out = new File('${tag}_genome.json')
    def genome = [:]
    genome.name = "${tag}"
    genome.meta = [:]
    if("${meta.shortName}" != "null") {
      genome.meta << ["shortName" : "${meta.shortName}"]
    }
    if("${meta.source}" != "null") {
      genome.meta << ["source" : "${meta.source}"]
    }
    if("${meta.release}" != "null") {
      genome.meta << ["release" : "${meta.release}"]
    }
    if("${meta.citation}" != "null") {
      genome.meta << ["citation" : "${meta.citation}"]
    }
    genome.blocks = []
    idx.eachLine { line ->
      if(line.toLowerCase() =~ /^(chr|[0-9]|x|y)/ ) {
        toks = line.split('\t')
        genome.blocks += [ "scope": toks[0].replaceFirst("^(C|c)(H|h)(R|r)[_]?",""), "featureType": "linear", "range": [1, toks[1].toInteger()] ]
      }
    }
    out.text = prettyPrint(toJson(genome))
    """
}


/*
* Given a FASTA with representative peps and the corresponding gtfgff3
* output FASTA with representative peps and definition lines
* mimicking the ensembl plants (EP) format for such data - this can then
* be piped into the same processes which we use for chewing through EP data
*/
process convertReprFasta2EnsemblPep {
  tag{tag}
  label 'fastx'

  input:
    //val arr from localInput
    set (val(meta), file(gtfgff3), file(reprPep)) from localInput

  output:
    set val(meta), file(pep) into localPepSeqs4Features, localPepSeqs4Aliases1, localPepSeqs4Aliases2

  script:
    tag=getAnnotationTagFromMeta(meta)
    //TRIAL RUN? ONLY TAKE FIRST n LINES
    cmd = trialLines != null ? "head -n ${trialLines}" : "cat"
    if(meta.containsKey("gtfgff3") && (gtfgff3.name).matches(".*gtf\$")) {
      // println("MATCHED gtf: "+gtfgff3)
        """
        ${cmd} ${reprPep} |  fasta_formatter | gtfAndRepr2ensembl_pep.awk -vversion="${meta.version}" - ${gtfgff3} > pep
        """
    } else if(meta.containsKey("gtfgff3") && (gtfgff3.name).matches(".*gff(3)?\$")) { //if(meta.containsKey("gff3")) {
      // println("MATCHED gff3: "+gtfgff3)
        """
        ${cmd} ${reprPep} | fasta_formatter | gff3AndRepr2ensembl_pep.awk -vversion="${meta.version}"  - ${gtfgff3} > pep
        """
    } else { //ASSUMING ENSEMBL PLANTS-LIKE FORMATTED PEPTIDE FASTA
      // println("NOT MATCHED gtfgff3: "+gtfgff3)
        """
        cp --no-dereference ${reprPep} pep
        """
    }
}

/*
 Only keep "representative" splice form for each gene,
 current approach selects longest transcript,
 we previously relied on ID suffix ".1", some times "-01"
 but some of the more recent Ensembl plants data sets
 no longer follow this convention
*/
process filterForRepresentativePeps {
  tag{meta}
  label 'fastx'
  input:
    set val(meta), file(pep) from remotePepSeqs

  output:
    set val(meta), file("${tag}_repr.pep") into remotePepSeqs4Features, remotePepSeqs4Aliases1, remotePepSeqs4Aliases2

  script:
    tag=getAnnotationTagFromMeta(meta)
    """
    fasta_formatter < ${pep} | paste - - | filterForRepresentative.awk > ${tag}_repr.pep
    [ -s ${tag}_repr.pep ] || (echo 'Error! Empty output file! ${tag}_repr.pep'; exit 1)
    """
}


/*
* Generate for pretzel JSON aliases linking features between chromosomes/genomes
* Fails if JSON invalid
*/
process generateFeaturesJSON {
  tag{tag}
  label 'json'
  label 'groovy'

  input:
    set val(meta), file(pep) from localPepSeqs4Features.mix(remotePepSeqs4Features)

  output:
    file "*.json.gz" into featuresJSON

  script:

    tag=getAnnotationTagFromMeta(meta)
    genome=getDatasetTagFromMeta(meta)
    shortName = (meta.containsKey("shortName") ? meta.shortName+"_genes" : "")
    shortName +=(meta.containsKey("annotation") ? "_"+meta.annotation : "") //only for cases where multiple annotations per genome
    """
    #!/usr/bin/env groovy

    import static groovy.json.JsonOutput.*
    pep = new File('${pep}').text
    out = new File('${tag}_annotation.json')
    def annotation = [:]
    annotation.meta = [:]
    if("${meta.shortName}" != "null") {
      annotation.meta << ["shortName" : "${shortName}"]
    }
    if("${meta.source}" != "null") {
      annotation.meta << ["source" : "${meta.source}"]
    }
    if("${meta.release}" != "null") {
      annotation.meta << ["release" : "${meta.release}"]
    }
    if("${meta.citation}" != "null") {
      annotation.meta << ["citation" : "${meta.citation}"]
    }
    annotation.name = "${tag}_genes"
    annotation.namespace = "${genome}:${tag}_annotation"
    annotation.parent = "${genome}"
    annotation.blocks = []
    TreeMap scope = [:] //keep keys sorted as the corresponding blocks get displayed in order in pretzel
    pep.eachLine { line ->
      if(line =~ /^>/ ) {
        toks = line.split()
        location = toks[2].split(":")
        gene = toks[3].split(":")
        key = location[2].replaceFirst("^(C|c)(H|h)(R|r)[_]?","")
        if(!scope.containsKey(key)) {
          scope << [(key) : []]
        }
        scope[key] << ["name" : gene[1], "range" : [ location[3].toInteger(), location[4].toInteger() ]]
      }
    }
    //GROUP TOGETHER FEATURES FROM/IN SAME BLOCK
    scope.each { k, features ->
      current = [ "scope": k, "featureType": "linear", "features": []]
      features.each { feature ->
        current.features << feature
      }
      annotation.blocks << current
    }
    out.text = prettyPrint(toJson(annotation))
    'gzip ${tag}_annotation.json'.execute()
    """
}


// //REPEAT INPUT FOR EACH SUBGENOME
// localPepSeqs4AliasesRep = Channel.create()
// localPepSeqs4Aliases.subscribe onNext: {
//   // println it[0]
//   if(it[0].containsKey("subgenomes")) {
//     for(subgenome in it[0].subgenomes) {
//       clone = it[0].clone()
//       clone.subgenome = subgenome
//       localPepSeqs4AliasesRep << [clone,it[1]]
//     }
//   } else {
//     localPepSeqs4AliasesRep << it
//   }
// }, onComplete: { localPepSeqs4Aliases.close(); localPepSeqs4AliasesRep.close() }


//COMBINE AND FILTER DUPLICATED CHANNEL TO ALLOW ALL VS ALL DATASETS COMPARISONS
// remotePepSeqs4AliasesCombined = remotePepSeqs4Aliases1.mix(localPepSeqs4AliasesRepSplit1).combine(remotePepSeqs4Aliases2.mix(localPepSeqs4AliasesRepSplit2))
pepSeqs4AliasesCombined = remotePepSeqs4Aliases1.mix(localPepSeqs4Aliases1).combine(remotePepSeqs4Aliases2.mix(localPepSeqs4Aliases2))
  .filter { getAnnotationTagFromMeta(it[0]) <= getAnnotationTagFromMeta(it[2])  }  //[species,version,file.pep]

// .collect().subscribe{ println it.combinations().each { a, b -> a[0].species < b[0].species} }

/*
* Identify best hit for each pep
*/
process pairProteins {
  tag{meta}
  label 'MMseqs2'
  errorStrategy 'ignore'

  input:
    set val(metaA), file('pepA'), val(metaB), file('pepB') from pepSeqs4AliasesCombined

  output:
     set val(metaA), val(metaB), file("*.tsv"), file(idlines) into pairedProteins

  script:
    tagA=getAnnotationTagFromMeta(metaA)
    tagB=getAnnotationTagFromMeta(metaB)
    meta = ["query": tagA, "target": tagB]
    basename=tagA+"_VS_"+tagB
    """
    mmseqs easy-search ${pepA} ${pepB} ${basename}.tsv \${TMPDIR:-/tmp}/${basename} \
    --format-mode 2 \
    -c ${params.minCoverage} \
    --min-seq-id ${params.minIdentity} \
    --threads ${task.cpus} -v 1 \
    && grep --no-filename '^>' ${pepA} ${pepB} | sed 's/^>//' > idlines
    """
    //'qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'
}

/*
* Generate JSON aliases linking features between chromosomes/genomes
*/
process generateAliasesJSON {
  tag{basename}
  label 'json'
  errorStrategy 'ignore'  //expecting exit status 3 if no aliases generated
  echo 'true'

  input:
    set(val(metaA), val(metaB), file(paired), file(idlines)) from pairedProteins

  output:
    file "${basename}_aliases.json.gz" into aliasesJSON
    // set val(basename), file("${basename}_aliases.json") into aliasesJSON

  script:
    // println(prettyPrint(toJson(metaA)))
    genome1=getDatasetTagFromMeta(metaA)
    genome2=getDatasetTagFromMeta(metaB)
    tag1=getAnnotationTagFromMeta(metaA)
    tag2=getAnnotationTagFromMeta(metaB)
    basename=tag1+"_VS_"+tag2
    namespace1=genome1+":"+tag1+"_annotation"
    namespace2=genome2+":"+tag2+"_annotation"
    cmd = tag1 != tag2 ? "cat ${paired} " : "excludeSameChromosome.awk -vtag1=${tag1} -vtag2=${tag2} ${idlines} ${paired}"
    """
    #at least one of the aligned pair must meet the minCoverageFilter threshold
    ${cmd} | awk '\$3 >= ${params.minIdentityFilter} && ((\$8-\$7+1)/\$13 >= ${params.minCoverageFilter} || (\$10-\$9+1)/\$14 >= ${params.minCoverageFilter})' \
    | blasttab2json.awk -vnamespace1=${namespace1} -vnamespace2=${namespace2} | gzip > ${basename}_aliases.json.gz
    zcat ${basename}_aliases.json | head | grep '[a-Z0-9]' > /dev/null || (echo "No aliases generated for ${basename}" && exit 3)
    """
}


process pack {
  label 'archive'
  executor 'local'

  input:
    file('*') from genomeBlocksJSON.collect()
    file('*') from featuresJSON.collect()
    file('*') from aliasesJSON.collect()

  output:
    file('*') into targzJSON

  """
  tar chzvf JSON-\$(date --iso-8601).tar.gz *.json *.json.gz
  tar chzvf  JSON-\$(date --iso-8601)-no_LC.tar.gz *.json *.json.gz --exclude '*LC*'
  """
}

// process mergeAliasesJSON {
//   label 'json'
//   tag{out}

//   input:
//     val tuple from aliasesJSON.groupTuple() //basename followed by a list of one or more *_alisaes.json file names

//   output:
//     file "*"

//   script:
//     out=tuple[0]
//     files=tuple[1].join(" ")
//     if(tuple[1].size() < 2) {
//       """
//       cp --preserve=links ${files} ${out}_aliases.json
//       """
//     } else {
//       """
//       echo "[" > ${out}_aliases.json
//       for f in ${files}; do
//         sed '1d;\$d' \${f} | sed 's/}\$/},/'
//       done | sed '\$d' >> ${out}_aliases.json
//       echo -e "    }\n]" >> ${out}_aliases.json
//       """
//     }
// }
