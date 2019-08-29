#!/usr/bin/env nextflow

if(!workflow.profile.contains('EP')) {
  println("This workflow requires -profile EP to be specified")
  exit 1
}

//INPUT PARAMS
trialLines = params.trialLines
// eprelease = params.eprelease

import static groovy.json.JsonOutput.*

//Otherwise JSON generation triggers stackoverflow when encountering Path objects
jsonGenerator = new groovy.json.JsonGenerator.Options()
                .addConverter(java.nio.file.Path) { java.nio.file.Path p, String key -> p.toUriString() }
                .build()


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
fastasuffix = params.fastasuffix

// def noKeyOrIsMap

//LOCAL marker sets

Channel.from(params.markers)
.map {
  [it, file(it.fasta)]
}
.set{ markerSetsChannel }

//LOCAL INPUTS
localInput = Channel.create()
localIndices = Channel.create()
localGenomeSeqs = Channel.create()
if(params.localAssembly != "NA") {
  params.localAssembly.each {
     //Genome Fasta optional (?)
    if(it.containsKey("fasta")) {
      localGenomeSeqs << [it,file(it.fasta)]
      it.remove('fasta') //preventing cached re runs after fasta added to meta
    }
    // println(prettyPrint(toJson(it)))
    // for (key in ["gtf","gff3"]) {
    key = 'gtfgff3'
    // if(it.containsKey(key)) {
    //IF MORE THAN ONE ANNOTATION PER GENOME
    if(it.containsKey("pep")) {
      if(it.pep instanceof Map) { // && it.containsKey(key) && it.get(key) instanceof Map) {
        it.pep.each {id, pep ->
        // println(id+" -> "+pep)
          clone = it.clone()
          clone.annotation = id
          gtfgff3 = (it.containsKey(key) && it.get(key).containsKey(id)) ? file(it.get(key).get(id)) : null
            localInput << [clone,gtfgff3,file(pep)]
        }
      } else { //if (!(it.pep instanceof Map) && !(it.get(key) instanceof Map)){
        gtfgff3 = it.containsKey(key) ? file(it.get(key)) : null
        localInput << [it,gtfgff3,file(it.pep)]
        // } else {
        //   exitOnInputMismatch(it)
      }
    }
    //ALL SHOULD HAVE AN INDEX
    if(it.containsKey("idx")) {
      localIndices << [it,file(it.idx)]
    }

  }
}
localInput.close()
localIndices.close()
localGenomeSeqs.close()

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
    set val(meta), file("${basename}.fasta") into remoteGenomeSeqs

  script:
    meta=["species":species, "version":version, "source": "https://plants.ensembl.org/"+species, "release": eprelease, "shortName": shortName]
    basename=getDatasetTagFromMeta(meta)
    idxurl=urlprefix+eprelease+"/fasta/"+species.toLowerCase()+"/dna_index/"+species+"."+version+idxsuffix
    fastaurl=urlprefix+eprelease+"/fasta/"+species.toLowerCase()+"/dna/"+species+"."+version+fastasuffix
    pepurl=urlprefix+eprelease+"/fasta/"+species.toLowerCase()+"/pep/"+species+"."+version+pepsuffix
    //Someone decided to embed Genus_species in version as well,
    //must be kept there to fetch from Ensembl plants but otherwise annoying as makes data set names long and repetitive
    //could be solved by explicitly setting pathis in input config (e.g. conf/triticeae.config)
    meta.version = ((meta.version-meta.species).strip('_'))
    if(trialLines == null) {
      """
      curl $idxurl > ${basename}.idx
      curl $fastaurl | gunzip --stdout > ${basename}.fasta
      curl $pepurl | gunzip --stdout > ${basename}.pep
      """
    } else {
      """
      curl $idxurl > ${basename}.idx
      curl $fastaurl | gunzip --stdout | head -n ${trialLines} > ${basename}.fasta
      curl $pepurl | gunzip --stdout | head -n ${trialLines} > ${basename}.pep
      """
    }
}

process alignMarkers {
  label 'minimap2'
  tag {"${refmeta} <- ${markersmeta}"}
  input:
    set val(refmeta), file(ref), val(markersmeta), file(markers) from remoteGenomeSeqs.mix(localGenomeSeqs).combine(markerSetsChannel).first() // <=========

  output:
    set val(outmeta), file('*.paf') into alignedMarkersChannel

  script:
  outmeta = [ref: refmeta, markers: markersmeta]
  """
  minimap2 -x sr \
  --secondary=yes \
  --cs=long \
  -I 30G -t ${task.cpus} ${ref} ${markers} > ${markers}_vs_${ref}.paf
  """
  /*
  Breakdown of -x sr Short single-end reads without splicing i.e.:
   -k21 Minimizer k-mer length [15]
   -w11 Minimizer window size [2/3 of k-mer length]. A minimizer is the smallest k-mer in a window of w consecutive k-mers.
   --sr Enable short-read alignment heuristics. In the short-read mode, minimap2 applies a second round of chaining with a higher minimizer occurrence threshold if no good chain is found. In addition, minimap2 attempts to patch gaps between seeds with ungapped alignment.

   -A2 Matching score [2]
   -B8 Mismatching penalty [4]
   -O12,32 Gap open penalty [4,24]. If INT2 is not specified, it is set to INT1.
   -E2,1 Gap extension penalty [2,1]. A gap of length k costs min{O1+k*E1,O2+k*E2}. In the splice mode, the second gap penalties are not used.
   -r50 Bandwidth used in chaining and DP-based alignment [500]. This option approximately controls the maximum gap size.
   -p.5 Minimal secondary-to-primary score ratio to output secondary mappings [0.8]. Between two chains overlaping over half of the shorter chain (controlled by -M), the chain with a lower score is secondary to the chain with a higher score. If the ratio of the scores is below FLOAT, the secondary chain will not be outputted or extended with DP alignment later. This option has no effect when -X is applied.
   -N20 Output at most INT secondary alignments [5]. This option has no effect when -X is applied.
   -f1000,5000 If fraction, ignore top FLOAT fraction of most frequent minimizers [0.0002]. If integer, ignore minimizers occuring more than INT1 times. INT2 is only effective in the --sr or -xsr mode, which sets the threshold for a second round of seeding.
   -n2 Discard chains consisting of <INT number of minimizers [3]
   -m20 Discard chains with chaining score <INT [40]. Chaining score equals the approximate number of matching bases minus a concave gap penalty. It is computed with dynamic programming.
   -s40 Minimal peak DP alignment score to output [40]. The peak score is computed from the final CIGAR. It is the score of the max scoring segment in the alignment and may be different from the total alignment score.
   -g200 Stop chain enlongation if there are no minimizers within INT-bp [10000].
   -2 Use two I/O threads during mapping. By default, minimap2 uses one I/O thread. When I/O is slow (e.g. piping to gzip, or reading from a slow pipe), the I/O thread may become the bottleneck. Apply this option to use one thread for input and another thread for output, at the cost of increased peak RAM.
   -K50m Number of bases loaded into memory to process in a mini-batch [500M]. Similar to option -I, K/M/G/k/m/g suffix is accepted. A large NUM helps load balancing in the multi-threading mode, at the cost of increased memory.
   --heap-sort=yes 	If yes, sort anchors with heap merge, instead of radix sort. Heap merge is faster for short reads, but slower for long reads. [no]
   --secondary=no Whether to output secondary alignments [yes]


   //IRRELEVANT FOR MARKERS AS CURRENTLY SET UP
   --frag=yes Whether to enable the fragment mode [no] ?????
   */
}

// alignedMarkersChannel.view { it -> groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it))}

process generateFeaturesFromMarkerAlignmentsJSON {
  tag { tag }
  label 'groovy'
  label 'json'

  input:
    set val(meta), file(paf) from alignedMarkersChannel

  output:
    file "*.json.gz" into markersJSON
    file "*.counts" into markersCounts

  script:
  tag=[meta.ref.species, meta.ref.version, meta.markers.name].join('_')
  genome=getDatasetTagFromMeta(meta.ref) //parent
  //shortName = (meta.ref.containsKey("shortName") ? meta.ref.shortName+"_"+meta.markers.name : "")

  // println(groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(meta)))

  template 'paf2pretzel.groovy'
  // """

  // """
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
    file "*.json" into genomeBlocksJSON, genomeBlocksStats

  script:
    tag=getDatasetTagFromMeta(meta)
    """
    #!/usr/bin/env groovy

    import static groovy.json.JsonOutput.*
    idx = new File('${idx}').text
    out = new File('${tag}_genome.json')
    def genome = [:]
    genome.name = "${tag}"
    genome.public = ${!params.makePrivate}
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
    genome.meta << ["type" : "Genome"]
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
 Given a FASTA with representative peps and the corresponding gtfgff3
 output FASTA with representative peps and definition lines
 mimicking the ensembl plants (EP) format for such data - this can then
 be piped into the same processes which we use for chewing through EP data
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
    file "*.counts" into featuresCounts

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
    counts = new File('${tag}_annotation.counts')
    def annotation = [:]
    annotation.public = ${!params.makePrivate}
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
        //Skip non-chromosome blocks
        if(key.toLowerCase() =~ /^(chr|[0-9]|x|y)/ ) {
          if(!scope.containsKey(key)) {
            scope << [(key) : []]
          }
          scope[key] << ["name" : gene[1], "value" : [ location[3].toInteger(), location[4].toInteger() ]]
        }
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
    //RECORD NUM FEATURS PER BLOCK
    counts.withWriterAppend{ wr ->
      annotation.blocks.each {
        wr.println  annotation.name+"\\t"+it.scope+"\\t"+it.features.size()
      }
    }
    //OUTPUT JSON, COMPRESS
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
  label 'jq'
  errorStrategy 'ignore'  //expecting exit status 3 if no aliases generated
  echo 'true'

  input:
    set(val(metaA), val(metaB), file(paired), file(idlines)) from pairedProteins

  output:
    file "${basename}_aliases.json.gz" into aliasesJSON
    file "${basename}_aliases.len" into aliasesCounts
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
    | blasttab2json.awk -vnamespace1=${namespace1} -vnamespace2=${namespace2} \
    | tee >(jq length >> ${basename}_aliases.len) \
    | gzip > ${basename}_aliases.json.gz \
    && zcat ${basename}_aliases.json | head | grep -E '[[:alnum:]]' > /dev/null || (echo "No aliases generated for ${basename}" && exit 3)
    """
}


process pack {
  label 'archive'
  executor 'local'

  input:
    file('*') from genomeBlocksJSON.collect()
    file('*') from featuresJSON.collect()
    file('*') from aliasesJSON.collect()
    file('*') from markersJSON.collect()

  output:
    file('*') into targzJSON

  """
  tar chzvf JSON-\$(date --iso-8601).tar.gz *.json *.json.gz
  tar chzvf  JSON-\$(date --iso-8601)-no_LC.tar.gz *.json *.json.gz --exclude '*LC*'
  """
}

process stats {
  echo true
  label 'summary'

  input:
    file('*') from genomeBlocksStats.collect()
    file('*') from featuresCounts.collect()
    file('*') from aliasesCounts.collect()
    file('*') from markersCounts.collect()

  output:
    file('*') into outstats

  script:
  """
  jq -r  '.blocks[] | (input_filename, .scope, .range[1])' *_genome.json | paste - - - | sort -V > blocks.counts
  cat *_annotation.counts | sort -V > feature.counts
  cat *_markers.counts | sort -V > markers.counts
  grep "" *_aliases.len > aliases.counts
  """
  //jq '.blocks[]' ${f} | jq 'input_filename, .scope, (.features | length)' | paste - - | sort -V
}

// process alignMarkers {
//   input:
//     file markers from markerSetsChannel

//   """
//   minimap2 -a -x sr -I 20G -t 20  161010_Chinese_Spring_v1.0_p
// seudomolecules.fasta 90kSNPprobes.seq.fasta > 90kSNPprobes.seq.fasta_vs_161010_Chinese_Spring_v1.0_pseudomolecules.sam
//   """
// }

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
