#!/usr/bin/env nextflow

if(workflow.profile.contains('BUSCOs')) {
  println("This workflow is not compatible with the BUSCOs profile")
  exit 1
}


// import static groovy.json.JsonOutput.*
//For pretty-printing nested maps etc
import groovy.json.JsonGenerator
import groovy.json.JsonSlurper
import groovy.json.JsonOutput

//Preventing stack overflow on Path objects and other  when map -> JSON
JsonGenerator jsonGenerator = new JsonGenerator.Options()
                .addConverter(java.nio.file.Path) { java.nio.file.Path p, String key -> p.toUriString() }
                .addConverter(Duration) { Duration d, String key -> d.durationInMillis }
                .addConverter(java.time.OffsetDateTime) { java.time.OffsetDateTime dt, String key -> dt.toString() }
                .addConverter(nextflow.NextflowMeta) { nextflow.NextflowMeta m, String key -> m.toJsonMap() }  //incompatible with Nextflow <= 19.04.0
                .excludeFieldsByType(java.lang.Class) // .excludeFieldsByName('class')
                // .excludeNulls()
                .build()

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


// def noKeyOrIsMap

//LOCAL marker/contigs to place sets

Channel.from(params.sequencesToPlace)
.map {
  [it, file(it.fasta)]
}
.set{ sequencesToPlaceChannel }

//INPUT DATA




/*
  Generic method for extracting a string tag or a file basename from a metadata map
 */
def getAnnotationTagFromMeta(meta, delim = '_') {
  return meta.species+delim+meta.version+(meta.containsKey("annotation") ? delim+meta.annotation : "")
}


/*
  Generic method for extracting a string tag or a file basename from a metadata map
 */
def getDatasetTagFromMeta(meta, delim = '_') {
  return meta.species+delim+meta.version
}

Channel.from(params.references)
  .into { refsChannel1; refsChannel2 ; refsChannel3}



process alignToGenome {
  label 'minimap2'
  tag {"${refmeta.subMap(['species','version'])} <- ${seqsmeta.name}"}

  input:
    tuple val(refmeta), file(ref), val(seqsmeta), file(seqs) from refsChannel2
        .filter { it.containsKey('fasta') }
        .map { [it, file(it.fasta)]}
        .combine(sequencesToPlaceChannel)

  output:
    set val(outmeta), file('*.paf') into alignedSeqsChannel

  when: //if target not specified place on all, otherwise check if current ref specified as target
    !seqsmeta.containsKey('target') || seqsmeta.target.any { it.species == refmeta.species && it.version == refmeta.version }


  script:
  outmeta = [ref: refmeta, seqs: seqsmeta.subMap(['name', 'seqtype'])]
  //preset: short read OR long assembly to ref OR long, HQ spliced to ref
  preset = seqsmeta.seqtype == 'markers' ? 'sr' : seqsmeta.seqtype == 'genomic' ? 'asm5' : 'splice:hq'
  secondary = seqsmeta.seqtype.toLowerCase().matches('markers|transcripts|cds|orf') ? 'yes' : 'no'
  csTag = seqsmeta.seqtype == 'markers' ? '--cs=long' : '' //save space by not printing cs in non-maraker modes
  alnParams = "-x ${preset} --secondary=${secondary} ${csTag} -I 30G"
  outmeta.align = [tool: 'minimap2', params: alnParams]
  """
  minimap2 ${alnParams} -t ${task.cpus} ${ref} ${seqs} > ${seqs}_vs_${ref}.paf
  """

  //TODO gzip/pigz --fast the output paf

  /*

  -N 50 \ ?????

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

// alignedSeqsChannel.view { prettyPrint(jsonGenerator.toJson(it)) }


process generateFeaturesFromSeqAlignmentsJSON {
  tag { tag }
  label 'groovy'
  label 'json'
  label 'tsv'
  label 'gff'
  label 'mem'
  errorStrategy { task.exitStatus = 3 ? 'ignore' : 'terminate' }  //expecting exit status 3 if no features placed which is valid e.g. when no good-enough alignments found

  input:
    set val(meta), file(paf) from alignedSeqsChannel

  output:
    file "*.{json.gz,tsv,gff}" optional true into placedSeqs
    file "*.counts" into placedSeqsCounts

  script:
  tag=[meta.ref.species, meta.ref.version, meta.seqs.name].join('_')
  genome=getDatasetTagFromMeta(meta.ref) //parent
  //shortName = (meta.ref.containsKey("shortName") ? meta.ref.shortName+"_"+meta.markers.name : "")

  // println(groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(meta)))

  // template 'paf2pretzel.groovy'
  """
  paf2pretzel.groovy \\
    --paf ${paf} \\
    --parent ${genome} \\
    --sequence-type ${meta.seqs.seqtype} \\
    --base-name ${tag} \\
    --namespace ${meta.seqs.name} \\
    --short-name ${meta.seqs.name} \\
    --align-tool ${meta.align.tool} \\
    --align-params "${meta.align.params}" \\
    --min-identity ${params.minMarkerIdentity} \\
    --min-coverage ${params.minMarkerCoverage} \\
    --allowed-target-id-pattern '${meta.ref.allowedIdPattern}' \\
    --output ${tag}_${meta.seqs.seqtype}.json.gz \\
    --out-tsv ${tag}_${meta.seqs.seqtype}.tsv \\
    --out-gff ${tag}_${meta.seqs.seqtype}.gff \\
    --out-counts ${tag}_${meta.seqs.seqtype}.counts
  """
}

refsChannel1
  .branch { meta -> //redirect data sets; ones without fai idx will need to have it generated
    ready: meta.containsKey('idx')
      [meta, file(meta.idx)]
    faidx: meta.containsKey('fasta')
      [meta, file(meta.fasta)]
  }
  .set { refs4genomeBlocks1 }

process faidxAssembly {
  tag{tag}
  label 'samtools'

  input:
    tuple val(meta), file(fasta) from refs4genomeBlocks1.faidx

  output:
    tuple val(meta), file("${fasta}.fai") into refs4genomeBlocks2

  script:
    tag=getDatasetTagFromMeta(meta)
    """
    #if err, likely due to gzipped not bgzipped fasta then index flat - we just need the lengths not the index!
    samtools faidx ${fasta} || (zcat ${fasta} > tmp && samtools faidx tmp && mv tmp.fai ${fasta}.fai)
    """
}

/*
* Generate genome blocks definitions JSON for pretzel
*/
process generateGenomeBlocksJSON {
  tag{tag}
  label 'json'
  label 'groovy'

  input:
    tuple val(meta), file(idx) from refs4genomeBlocks1.ready.mix(refs4genomeBlocks2)

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
      if(line.toLowerCase() =~ /^(ch|[0-9]{1,2}|x|y|i|v)/ || line ==~ '${meta.allowedIdPattern}' ) {
        toks = line.split('\\t| ')
        genome.blocks += [ "scope": toks[0].replaceFirst("^(C|c)(H|h)(R|r)[_]?",""), "featureType": "linear", "range": [1, toks[1].toInteger()] ]
      }
    }
    if(genome.blocks.isEmpty()) {
      System.err.println('No blocks defined for ${tag}, this may be caused by chromosome naming, terminating')
      System.exit(2)
    }
    out.text = prettyPrint(toJson(genome))
    """
}

refsChannel3
  // .view {
  //   """
  //   ${it}
  //   ${it.containsKey('gff3')}
  //   ${it.containsKey('gtf')}
  //   ${(it.containsKey('gff3') || it.containsKey('gtf'))}
  //   ${!(it.containsKey('gff3') || it.containsKey('gtf'))}
  //   """
  // }
  .filter { meta -> meta.containsKey('pep') }
  .map { meta -> // DUPLICATE EMISSIONS IF MULTIPLE ANNOTATIONS PER REFERENCE ASSEMBLY
    if(meta.pep instanceof Map) {
      def repeated = []
      meta.pep.each { k,v ->
        def item = meta.subMap(meta.keySet().minus(['pep','gff3','gtf'])) + [pep: v, annotation: k]
        if(meta.containsKey('gff3') && meta.gff3.containsKey(k)) {
          item.gff3 = meta.gff3."${k}"
        } else if(meta.containsKey('gtf') && meta.gtf.containsKey(k)) {
          item.gtf = meta.gtf."${k}"
        }
        repeated << item
      }
      repeated
    } else {
      meta
    }
  }
  .flatten()
  .branch { meta -> //redirect data sets; ones with pep but without gff/gtf are assumed to be in ENSEMBL format
     pep4Conversion: (meta.containsKey('gff3') || meta.containsKey('gtf'))
     pepEnsembl: !(meta.containsKey('gff3') && !meta.containsKey('gtf'))
       [meta, file(meta.pep)]
  }
  .set { refsWithPepChannel }

// refsWithPepChannel.pep4Conversion.view { it -> groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it))}

/*
 Only keep "representative" splice form for each gene,
 current approach selects longest transcript,
 we previously relied on ID suffix ".1", some times "-01"
 but some of the more recent Ensembl plants data sets
 no longer follow this convention
*/
process filterForRepresentativePeps {
  tag{meta.subMap(['species','version'])}
  label 'fastx'
  input:
    set val(meta), file(pep) from refsWithPepChannel.pepEnsembl

  output:
    set val(meta), file("${tag}_repr.pep.gz") into representativePepSeqs4Features, representativePepSeqs4Aliases1, representativePepSeqs4Aliases2

  script:
    tag=getAnnotationTagFromMeta(meta)
    cmd = "${pep}".endsWith(".gz") ? "zcat" : "cat"
    """
    ${cmd} ${pep} | fasta_formatter | paste - - | filterForRepresentative.awk | gzip -c > ${tag}_repr.pep.gz
    [ ! -z \$(zcat ${tag}_repr.pep.gz | head -c1) ] || (echo 'Error! Empty output file! ${tag}_repr.pep.gz'; exit 1)
    """
}


/*
 Given a FASTA with representative peps and the corresponding gtfgff3
 output FASTA with representative peps and definition lines
 mimicking the ensembl plants (EP) format for such data - this can then
 be piped into the same processes which we use for chewing through EP data
*/
process convertReprFasta2EnsemblPep { //TODO - NOT WORKING IF ENSEMB-FORMATTED INPUT (should not be used here but need to pass-through if already formatted?)
  tag{tag}
  label 'fastx'

  input:
    tuple (val(meta), file(gtfgff3), file(reprPep)) from refsWithPepChannel.pep4Conversion
      //.filter { meta -> meta.containsKey('pep') && (meta.containsKey('gff3') || meta.containsKey('gtf'))}
      .map { meta ->  [ meta, file( meta.containsKey('gff3') ? meta.gff3 : meta.gtf ), file( meta.pep ) ] }

  output:
    tuple val(meta), file('pep.gz') into pepSeqs4Features, pepSeqs4Aliases1, pepSeqs4Aliases2

  script:
    tag=getAnnotationTagFromMeta(meta)
    //TRIAL RUN? ONLY TAKE FIRST n LINES


    cmd0 = "${reprPep}".endsWith(".gz") ? "zcat" : "cat"
    cmd1 = "${gtfgff3}".endsWith(".gz") ? "zcat" : "cat"
    // if(meta.containsKey("gtfgff3") && (gtfgff3.name).matches(".*gtf\$")) {
    if(meta.containsKey("gtf")) {
        """
        ${cmd0} ${reprPep} |  fasta_formatter | gtfAndRepr2ensembl_pep.awk -vversion="${meta.version}" - <(${cmd1} ${gtfgff3}) | gzip  > pep.gz
        [ ! -z \$(zcat pep.gz | head -c1) ] || (echo 'Error! Empty output file! pep.gz'; exit 1)
        """
    } else { //if(meta.containsKey("gtfgff3") && (gtfgff3.name).matches(".*gff(3)?\$")) { //if(meta.containsKey("gff3")) {
      // println("MATCHED gff3: "+gtfgff3)
        """
        ${cmd0} ${reprPep} | fasta_formatter | gff3AndRepr2ensembl_pep.awk -vversion="${meta.version}"  - <(${cmd1} ${gtfgff3}) | gzip > pep.gz
        [ ! -z \$(zcat pep.gz | head -c1) ] || (echo 'Error! Empty output file! pep.gz'; exit 1)
        """
    // } else { //ASSUMING ENSEMBL PLANTS-LIKE FORMATTED PEPTIDE FASTA
    //   // println("NOT MATCHED gtfgff3: "+gtfgff3)
    //     """
    //     cp --no-dereference ${reprPep} pep
    //     """
    }
}




/*
* Generate for pretzel JSON aliases linking features between chromosomes/genomes
* Fails if JSON invalid
*/
process generateFeaturesJSON {
  tag{tag}
  label 'json'
  label 'groovy'
  echo true
  errorStrategy 'terminate'

  input:
    set val(meta), file(pep) from representativePepSeqs4Features.mix(pepSeqs4Features)
    // set val(meta), file(pep) from refsChannel2.map { meta -> [meta, file(meta.pep)] }

  output:
    file "*.json.gz" into featuresJSON
    file "*.counts" into featuresCounts

  script:
    tag=getAnnotationTagFromMeta(meta)
    genome=getDatasetTagFromMeta(meta)
    shortName = (meta.containsKey("shortName") ? meta.shortName+"_genes" : "")
    shortName +=(meta.containsKey("annotation") ? "_"+meta.annotation : "") //only for cases where multiple annotations per genome
    // """
    // ls -la
    // """
    """
    #!/usr/bin/env groovy

    import java.util.zip.GZIPInputStream
    import java.util.zip.GZIPOutputStream
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
    def pepStream = new FileInputStream(new File('${pep}'))
    def inStream = '${pep}'.endsWith('.gz') ? new GZIPInputStream(pepStream , 1024) : pepStream
    def content = new BufferedReader(new InputStreamReader(inStream, "UTF-8"), 1024);
    while ((line = content.readLine()) != null && !line.isEmpty() ) {
    // pep.eachLine { line ->
      if(line =~ /^>/ ) {
        toks = line.split()
        location = toks[2].split(":")
        gene = toks[3].split(":")
        key = location[2].replaceFirst("^(C|c)(H|h)(R|r)?[_]?","")
        //Skip non-chromosome blocks
        if(key.toLowerCase() =~ /^(ch|[0-9]|x|y|i|v)/ || key ==~ '${meta.allowedIdPattern}' ) {
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
    'gzip ${tag}_annotation.json'.execute().waitFor()
    """
}

pepSeqs4AliasesCombined = representativePepSeqs4Aliases1.mix(pepSeqs4Aliases1).combine(representativePepSeqs4Aliases2.mix(pepSeqs4Aliases2))
  .filter { getAnnotationTagFromMeta(it[0]) <= getAnnotationTagFromMeta(it[2])  }  //[species,version,file.pep]
  // .first()
  // .view{ [it[0].species, it[2].species] }
  // .view { it -> groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it))}

/*
* Identify best hit for each pep
*/
process pairProteins {
  tag{meta}
  label 'MMseqs2'
  errorStrategy 'ignore'

  input:
    set val(metaA), file('pepA.gz'), val(metaB), file('pepB.gz') from pepSeqs4AliasesCombined

  output:
     set val(metaA), val(metaB), file("*.tsv"), file(idlines) into pairedProteins

  script:
    tagA=getAnnotationTagFromMeta(metaA)
    tagB=getAnnotationTagFromMeta(metaB)
    meta = ["query": tagA, "target": tagB]
    basename=tagA+"_VS_"+tagB
    """
    mmseqs easy-search pepA.gz pepB.gz ${basename}.tsv \${TMPDIR:-/tmp}/${basename} \
    --format-mode 2 \
    -c ${params.minCoverage} \
    --min-seq-id ${params.minIdentity} \
    --threads ${task.cpus} -v 1 \
    && zcat pepA.gz pepB.gz | grep --no-filename '^>'  | sed 's/^>//' > idlines
    """
    //'qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'
}

/*
* Generate JSON aliases linking features between chromosomes/genomes
*/
process generateAliasesJSON {
  // maxForks 1
  tag{basename}
  label 'json'
  label 'jq'
  errorStrategy { task.exitStatus = 3 ? 'ignore' : 'terminate' }  //expecting exit status 3 if no aliases generated which is valid e.g. when dataset consisting of a single chr/block aligned to itself
  // errorStrategy 'finish'
  // echo 'true'

  input:
    set(val(metaA), val(metaB), file(paired), file(idlines)) from pairedProteins

  output:
    file "${basename}*_aliases.json.gz" optional true into aliasesJSON
    file "${basename}_aliases.len" into aliasesCounts

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
    # No aliases if all genes on a single chromososme, same assembly
    if [ \$(cut -f2,3 -d':' ${idlines} | sort |uniq | wc -l) -eq 1 ]; then
      echo 0 > ${basename}_aliases.len
      exit 0
    fi
    #at least one of the aligned pair must meet the minCoverageFilter threshold
    #and we split in to max 500k chunks to limit output JSON size
    ${cmd} | awk '\$3 >= ${params.minIdentityFilter} && ((\$8-\$7+1)/\$13 >= ${params.minCoverageFilter} || (\$10-\$9+1)/\$14 >= ${params.minCoverageFilter})' \
    | split -d -l 500000 - __part \
    && for PART in __part*; do
      blasttab2json.awk -vnamespace1=${namespace1} -vnamespace2=${namespace2} \${PART} \
      | tee >(jq length >> ${basename}_\${PART#__part}_aliases.len) \
      | gzip > ${basename}_\${PART#__part}_aliases.json.gz
    done \
    && zcat ${basename}*aliases.json.gz | head | grep -E '[[:alnum:]]' > /dev/null || (echo "No aliases generated for ${basename}" && exit 3) \
    && awk '{tot+=\$1};END{print tot}' ${basename}_*_aliases.len > ${basename}_aliases.len \
    && if [ \$(ls -1 __part*  | wc -l) -eq 1 ]; then mv ${basename}_\${PART#__part}_aliases.json.gz ${basename}_aliases.json.gz; fi \
    && rm __part*
    """
        //tmp: | split -dl5 --additional-suffix '.json' - part
    //    sed -e '1 i\[' -e '$ i\]'
}

process stats {
  label 'summary'
  label 'jq'

  input:
    // file('*') from genomeBlocksStats.collect()
    //                  .mix(featuresCounts.collect())
    //                  .mix(aliasesCounts.collect())
    //                  .mix(placedSeqsCounts.collect())
    file('*') from genomeBlocksStats.collect()
    file('*') from featuresCounts.collect()
    file('*') from aliasesCounts.collect()
    file('*') from placedSeqsCounts.collect()

  output:
    file('*') into outstats

  script:
  """
  jq -r  '.blocks[] | (input_filename, .scope, .range[1])' *_genome.json | paste - - - | sort -V > blocks.counts
  cat *_annotation.counts | sort -V > feature.counts
  cat *_{markers,transcripts,cds,genomic}.counts | sort -V > placed.counts
  grep "" *_aliases.len > aliases.counts
  echo "one or more of the above can fail - only relevant if not expected to..."
  """
  //jq '.blocks[]' ${f} | jq 'input_filename, .scope, (.features | length)' | paste - - | sort -V
}

process pack {
  label 'archive'
  executor 'local'

  input:
    file('*') from genomeBlocksJSON.collect()
    file('*') from featuresJSON.collect()
    file('*') from aliasesJSON.collect()
    file('*') from placedSeqs.collect()
    file('*') from outstats

  output:
    file('*') into targzJSON

  """
  tar chzvf JSON-\$(date --iso-8601).tar.gz *.json *.json.gz *.counts *.tsv *.gff
  tar chzvf  JSON-\$(date --iso-8601)-no_LC.tar.gz *.json *.json.gz *.counts *.tsv --exclude '*LC*'
  """
}

