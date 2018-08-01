#!/usr/bin/env nextflow
// echo true

//INPUT PARAMS
trialLines = params.trialLines
eprelease = params.eprelease

//LIST ELEMS ASSUMED TO MATCH
assemblySpeciesList = params.assemblySpecies.tokenize(",")
assemblyVersionList = params.assemblyVersion.tokenize(",")

//STATIC(?)
urlprefix = params.urlprefix
pepsuffix = params.pepsuffix
idxsuffix = params.idxsuffix

//LOCAL Zavitan for now...
// Channel.from(params.localAssemblySpecies.tokenize(","))
// .merge(Channel.from(params.localAssemblyVersion.tokenize(",")))
// .into { localIds1; localIds2 }

// localIds1.merge(Channel.fromPath(params.localgtfgff3)).merge(Channel.fromPath(params.localPep)).set{ localInputGtfGff3Pep }
// localIds2.merge(Channel.fromPath(params.localIdx)).set { localIndices }

//ARRANGE INPUTS FOR PROCESSES
localInputGtfGff3Pep = Channel.create()
localIndices = Channel.create()
if(params.localAssembly != "NA") {
  params.localAssembly.each {
    //EXPECT TO HAVE SOME DATASETS WITH gtfgff3, other with gff3 instead
    if(it.containsKey("gtf")) {
      localInputGtfGff3Pep << [it,file(it.gtf),file(it.pep)]
    } else if(it.containsKey("gff3")) {
      localInputGtfGff3Pep << [it,file(it.gff3),file(it.pep)]
    }
    //ALL SHOULD HAVE AN INDEX
    if(it.containsKey("idx")) {
      localIndices << [it,file(it.idx)]
    }
  }
}
localInputGtfGff3Pep.close()
localIndices.close()


// //ARRANGE INPUTS FOR PROCESSES
// referencesLocal = Channel.create()
// referencesRemote = Channel.create()
// params.references.each {
//   //EXPECT TO HAVE SOME DATASETS WITH fasta  
//   if(it.containsKey("fasta")) {
//     if((it.fasta).matches("^(https?|ftp)://.*\$")) {
//       referencesRemote << it
//     } else {
//       referencesLocal << [it,file(it.fasta)]
//     }
//   }
// }
// referencesRemote.close()
// referencesLocal.close()



/* 
  Generic method for extracting a string tag or a file basename from a metadata map
 */
def getTagFromMeta(meta, delim = '_') {
  return meta.species+delim+meta.version+(trialLines == null ? "" : delim+trialLines+delim+"trialLines")
}

/* 
  Generic method for extracting a string tag or a file basename from a metadata map allowing for an optional subGenome suffix
 */
def getUniqId(meta, delim = '_') { 
  return getTagFromMeta(meta, delim)+(meta.containsKey("subgenome") ? delim+meta.subgenome : "")  }



/*
* Download peptide seqs and assembly index files from Ensembl plants
*/
process fetchRemoteDataFromEnsemblPlants {
  tag{meta}
  label 'download'

  input:
    val species from assemblySpeciesList
    val version from assemblyVersionList
    // val urlprefix
    // val pepsuffix
    // val idxsuffix

  output:
    set val(meta), file("${basename}.idx") into remoteIndices
    set val(meta), file("${basename}.pep") into remotePepSeqs

  script:
    meta=["species":species, "version":version, "source": "https://plants.ensembl.org/"+species, "release": params.eprelease]
    basename=getTagFromMeta(meta)
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


// process fetchRemoteData {
//   tag{meta.subMap(['species','version'])}
//   label 'download'
//   echo true

//   input:
//     val(meta) from inputRemote

//   // output:
//   //   set val(meta), file("${basename}.fasta") into referencesRemoteFasta    

//   script:
//     basename=getTagFromMeta(meta)
//     cmd = ''
//     for(key in ['idx','pep','gtf','gff3']) {      
//       if(meta.containsKey(key)) {
//         value = meta.get(key)
//         println(key+" "+value)
//       }
//     }
//     """
//     echo
//     """
//     // //DECOMPRESS?
//     // cmd = ""
//     // cmd += (meta.fasta).matches("^.*\\.gz\$") ?  "| gunzip --stdout " :  " "
//     // cmd += (meta.fasta).matches("^.*\\.zip\$") ?  "| funzip " :  " "
//     // //TRIAL RUN? ONLY TAKE FIRST n LINES
//     // cmd += trialLines != null ? "| head -n ${trialLines}" : ""    
//     // """
//     // curl ${meta.fasta} ${cmd} > ${basename}.fasta
//     // """
// }


/*
* Generate genome blocks definitions JSON for pretzel
* Fails if JSON invalid
*/
process generateGenomeBlocksJSON {
  tag{tag}
  label 'json'
  label 'groovy'

  echo true

  input:
    set val(meta), file(idx) from localIndices.mix(remoteIndices)

  output:
    file "*.json" into genomeBlocksJSON

  script:
    tag=getTagFromMeta(meta)
    """
    #!/usr/bin/env groovy

    import static groovy.json.JsonOutput.*
    idx = new File('${idx}').text
    out = new File('${tag}_genome.json')
    def genome = [:] 
    genome.name = "${tag}"
    genome.meta = [:]
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
      if(line.toLowerCase() =~ /^(chr|[0-9])/ ) {
        toks = line.split('\t')
        genome.blocks += [ "scope": toks[0].replaceFirst("^(C|c)(H|h)(R|r)[_]?",""), "featureType": "linear", "range": [1, toks[1].toInteger()] ]
      }      
    }
    out.text = prettyPrint(toJson(genome))
    """

    // """
    // awk 'BEGIN{IGNORECASE=1;} \$1 ~/^(chr|[0-9])/' "${idx}" \
    // | faidx2json.awk -vname=${tag} \
    // | python -mjson.tool > ${tag}_genome.json
    // """
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
    set (val(meta), file(gtfgff3), file(reprPep)) from localInputGtfGff3Pep

  output:
    set val(meta), file(pep) into localPepSeqs4Features, localPepSeqs4Aliases1, localPepSeqs4Aliases2
  script:
    tag=getTagFromMeta(meta)    
    //TRIAL RUN? ONLY TAKE FIRST n LINES
    cmd = trialLines != null ? "head -n ${trialLines}" : "cat" 
    //cmd = "cat"
    if(meta.containsKey("gtf")) {
      """
      ${cmd} ${reprPep} |  fasta_formatter | gtfAndRepr2ensembl_pep.awk -vversion="${meta.version}" - ${gtfgff3} > pep
      """
    } else if(meta.containsKey("gff3")) {
      """
      ${cmd} ${reprPep} | fasta_formatter | gff3AndRepr2ensembl_pep.awk -vversion="${meta.version}"  - ${gtfgff3} > pep
      """
    }
}

/*
* Only keep "representative" splice form for each gene, typically suffix ".1", some times "-01"
*/
process filterForRepresentativePeps {
  tag{meta}
  //label 'fastx'
  input:
    set val(meta), file(pep) from remotePepSeqs

  output:
    set val(meta), file("${tag}_repr.pep") into remotePepSeqs4Features, remotePepSeqs4Aliases1, remotePepSeqs4Aliases2

  script:
  tag=getTagFromMeta(meta)

    """
    cat ${pep} | filterForRepresentative.awk > ${tag}_repr.pep
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
    //WARNINIG! HC OUTPUT CURRENTLY OVERWRITEN BY LC OUTPUT (ORE VICE-VERSA)

  output:
    file "*.json" into featuresJSON

  script:
    tag=getTagFromMeta(meta)
    // base=map.subgenomes.size() < 2 ? tag : tag+"_"+subgenome
    """
    #!/usr/bin/env groovy

    import static groovy.json.JsonOutput.*
    pep = new File('${pep}').text
    out = new File('${tag}_annotation.json')
    def annotation = [:] 
    annotation.meta = [:]
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
    annotation.namespace = "${tag}:${tag}_annotation"
    annotation.parent = "${tag}"
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
    """
    // awk '\$1 ~ /^>/' ${pep}  | sort -k3,3V | pep2featuresJSON.awk -vname="${tag}" \
    // | python -mjson.tool > ${tag}_annotation.json
    // """
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


// localPepSeqs4AliasesFeedback = Channel.create()
// localPepSeqs4AliasesFeedbackLoop = localPepSeqs4AliasesFeedback.filter { it[0].containsKey("subgenome") }

// process repeatPolyPeps{
//   echo true 

//   input:  
//     set val(meta), file(pep) from localPepSeqs4Aliases.mix(localPepSeqs4AliasesFeedbackLoop)
//   output:
//     set val(clonemeta), file(outpep) into localPepSeqs4AliasesRep
//     set val(repmeta), file(repep) into localPepSeqs4AliasesFeedback

//   // when:
//   //   !meta.containsKey("final")  

//   script:
//     // if(meta.containsKey("subgenomes")) {  
//     //   i = meta.containsKey("subgenome") ? meta.subgenomes.findIndexOf{it == meta.subgenome} : 0    
//     //   if(++i < meta.subgenomes.size()) {        
//     //     println (i+" "+meta.subgenomes[i])
//     //     clonemeta = meta.clone()
//     //     clonemeta.subgenome = meta.subgenomes[i]
//     //     repmeta = clonemeta.clone()
//     //     outpep = pep
//     //     repep = pep      
//     //   }
//     // } else {
//       clonemeta = meta
//       outpep = pep
//       repmeta = meta    
//       repmeta.final = true
//       repep = pep
//     // }
//     """
//     ls -lth
//     """
//   // }
// }

// Y = Channel.create()
// X = Channel.from([1,null],[2,["A","B"]],[3,["A","B","C"]]).map { it -> it }
//     .subscribe { println it };
// X.subscribe onNext: { 
//   println "BEFORE: "+it 
//     if(it[1]==null) {
//       Y << it
//     } else {
//       for(i in it[1]) {
//         clone = it.clone()
//         clone[1]=i
//         Y << clone
//       }
//     }
// }, onComplete: { X.close(); Y.close() }
// Y.subscribe { println "AFTER: "+ it}
//  Channel
//    .from([1,[species:"a"]],[2,[species: "b", sub:["A","B"]]],[3,[species:"c", sub:["A","B","C"]]])
//    .println()


// process splitPepSeqsPerSubGenome {
//   tag{tag}
//   echo true

//   input:
//     set val(meta), file(pep) from localPepSeqs4AliasesRep
    
//   output: //optional true to account for empty files produced during trial runs, e.g. no peps for subgenomes B, D but A only
//     set val(meta), file("${tag}.pep") into localPepSeqs4AliasesRepSplit1, localPepSeqs4AliasesRepSplit2 
//   //   // set val(map), file(pep) optional true into localPepSeqs4AliasesNoSubgenomes

//   script:    
//     tag=getUniqId(meta) //getTagFromMeta(meta)+(meta.containsKey("subgenome") ? "_"+meta.subgenome : "")
//     if(meta.containsKey("subgenome") && meta.subgenomes.size() > 1)  {
//        """
//        awk '{
//         if(\$1 ~ /^>/) {
//           split(\$3,arr,":");
//           if(arr[3] ~ /${meta.subgenome}\$/) {
//             current=1; print
//           } else {
//             current=0
//           }
//         } else if(current) {
//           print
//         }
//       }' pep > ${tag}.pep
//       """
//     } else {
//       """
//       cp --preserve=links pep ${tag}.pep
//       """
//     }
  
//   //   //>TRIDC1AG000070.2 pep chromosome:WEWseq_PGSB_20160501:chr1A:409955:411146 gene:TRIDC1AG000070

// }



//COMBINE AND FILTER DUPLICATED CHANNEL TO ALLOW ALL VS ALL DATASETS COMPARISONS
// remotePepSeqs4AliasesCombined = remotePepSeqs4Aliases1.mix(localPepSeqs4AliasesRepSplit1).combine(remotePepSeqs4Aliases2.mix(localPepSeqs4AliasesRepSplit2))
pepSeqs4AliasesCombined = remotePepSeqs4Aliases1.mix(localPepSeqs4Aliases1).combine(remotePepSeqs4Aliases2.mix(localPepSeqs4Aliases2))
  .filter { getUniqId(it[0]) <= getUniqId(it[2])  }  //[species,version,file.pep]

  // remotePepSeqs4AliasesCombined.subscribe {println (getUniqId(it[0])+" "+getUniqId(it[2])) }


// .collect().subscribe{ println it.combinations().each { a, b -> a[0].species < b[0].species} }


// [
  // [species:Arabidopsis_thaliana, version:TAIR10],
  // /home/rad/repos/pretzel-input-generator/work/59/e1b177dcc513678af115b00116e9b6/Arabidopsis_thaliana_TAIR10.pep,
  // [species:Brachypodium_distachyon, version:v1.0],
  // /home/rad/repos/pretzel-input-generator/work/82/761d637387c59c8b6ceeb50ea870e2/Brachypodium_distachyon_v1.0.pep]

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

  // when:
  //   !pepA.isEmpty() && !pepB.isEmpty()

  script:    
    // x = file(pepA)
    // println(x.getName())    
    // tagA=metaA.species+"_"+metaA.version //no subgenome spec to be pass to next process as used directly in JSON
    // tagB=metaB.species+"_"+metaB.version     
    // tag=tagA+(metaA.containsKey("subgenome") ? "_"+metaA.subgenome : "")+"_VS_"+tagB+(metaB.containsKey("subgenome") ? "_"+metaB.subgenome : "")
    // tagA=getTagFromMeta(metaA)
    // tagB=getTagFromMeta(metaB)
    meta = ["query": getTagFromMeta(metaA), "target": getTagFromMeta(metaB)]
    basename=getUniqId(metaA)+"_VS_"+getUniqId(metaB)
    // labtag=metaA.toString()+" VS "+metaB.toString()
    //TODO: EXPOSE PARAMS?    

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

// //  remotePepSeqs4Aliases1.combine(remotePepSeqs4Aliases2).filter { it[0] != it [3]  && it[0]+it[1] < it[3]+it[4]} .subscribe { println  "NOPE: $it" }



// import static groovy.json.JsonOutput.*  

/*
* Generate JSON aliases linking features between chromosomes/genomes
*/
process generateAliasesJSON {
  tag{basename}
  label 'json'
  errorStrategy 'ignore'

  input:
    set(val(metaA), val(metaB), file(paired), file(idlines)) from pairedProteins

  output:
    set val(outtag), file("${basename}_aliases.json") into aliasesJSON

  script:  
    // println(prettyPrint(toJson(metaA)))
    tag1=getTagFromMeta(metaA)
    tag2=getTagFromMeta(metaB)
    basename=getUniqId(metaA)+"_VS_"+getUniqId(metaB)
    outtag=tag1+"_VS_"+tag2
    namespace1=tag1+":"+tag1+"_annotation"
    namespace2=tag2+":"+tag2+"_annotation"    
    cmd = tag1 != tag2 ? "cat ${paired} " : "excludeSameChromosome.awk -vtag1=${tag1} -vtag2=${tag2} ${idlines} ${paired}"
    """    
    #at least one of the aligned pair must meet the minCoverageFilter threshold
    ${cmd} | awk '\$3 >= ${params.minIdentityFilter} && ((\$8-\$7+1)/\$13 >= ${params.minCoverageFilter} || (\$10-\$9+1)/\$14 >= ${params.minCoverageFilter})' \
    | blasttab2json.awk -vnamespace1=${namespace1} -vnamespace2=${namespace2} > ${basename}_aliases.json
    # python -mjson.tool > ${basename}_aliases.json
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
