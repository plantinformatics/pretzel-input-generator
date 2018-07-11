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
* Download peptide seqs and assembly index files from Ensembl plants
*/
process fetchRemoteDataFromEnsemblPlants {
  tag{map}
  label 'download'

  input:
    val species from assemblySpeciesList
    val version from assemblyVersionList
    // val urlprefix
    // val pepsuffix
    // val idxsuffix

  output:
    set val(map), file(idx) into remoteIndices
    set val(map), file(pep) into remotePepSeqs

  script:
    map=["species":species, "version":version]
    idxurl=urlprefix+eprelease+"/fasta/"+species.toLowerCase()+"/dna_index/"+species+"."+version+idxsuffix
    pepurl=urlprefix+eprelease+"/fasta/"+species.toLowerCase()+"/pep/"+species+"."+version+pepsuffix
    if(trialLines == null) {
      """
      curl $idxurl > idx
      curl $pepurl | gunzip --stdout > pep
      """
    } else {
      """
      curl $idxurl > idx
      curl $pepurl | gunzip --stdout | head -n ${trialLines} > pep
      """
    }
}


process fetchRemoteData {
  tag{meta.subMap(['species','version'])}
  label 'download'
  echo true

  input:
    val(meta) from inputRemote

  // output:
  //   set val(meta), file("${basename}.fasta") into referencesRemoteFasta    

  script:
    basename=getTagFromMeta(meta)
    cmd = ''
    for(key in ['idx','pep','gtf','gff3']) {      
      if(meta.containsKey(key)) {
        value = meta.get(key)
        println(key+" "+value)
      }
    }
    """
    echo
    """
    // //DECOMPRESS?
    // cmd = ""
    // cmd += (meta.fasta).matches("^.*\\.gz\$") ?  "| gunzip --stdout " :  " "
    // cmd += (meta.fasta).matches("^.*\\.zip\$") ?  "| funzip " :  " "
    // //TRIAL RUN? ONLY TAKE FIRST n LINES
    // cmd += trialLines != null ? "| head -n ${trialLines}" : ""    
    // """
    // curl ${meta.fasta} ${cmd} > ${basename}.fasta
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
    set (val(map), file(gtfgff3), file(reprPep)) from localInputGtfGff3Pep

  output:
    set val(map), file(pep) into localPepSeqs4Features, localPepSeqs4Aliases
  script:
    tag=map.species+"_"+map.version
    if(map.containsKey("gtf")) {
      """
      fasta_formatter < ${reprPep} | gtfAndRepr2ensembl_pep.awk -vversion="${map.version}" - ${gtfgff3} > pep
      """
    } else if(map.containsKey("gff3")) {
      """
      fasta_formatter < ${reprPep} | gff3AndRepr2ensembl_pep.awk -vversion="${map.version}"  - ${gtfgff3} > pep
      """
    }
}

/*
* Generate genome blocks definitions JSON for pretzel
* Fails if JSON invalid
*/
process generateGenomeBlocksJSON {
  tag{tag}
  label 'json'

  input:
    set val(map), file(idx) from localIndices.mix(remoteIndices)

  output:
     file "*.json"

  script:
    tag=map.species+"_"+map.version
    """
    awk '\$1 ~/^(chr|[0-9])/' "${idx}" \
    | faidx2json.awk -vname="${tag}" \
    | python -mjson.tool > "${tag}"_genome.json
    """
}

/*
* Only keep "representative" splice form for each gene, typically suffix ".1", some times "-01"
*/
process filterForRepresentativePeps {
  tag{map}
  //label 'fastx'
  input:
    set val(map), file(pep) from remotePepSeqs

  output:
    set val(map), file("${tag}.pep") into remotePepSeqs4Features, remotePepSeqs4Aliases1, remotePepSeqs4Aliases2

  script:
  tag=map.species+"_"+map.version

    """
    cat ${pep} | filterForRepresentative.awk > "${tag}.pep"
    """
}


/*
* Generate for pretzel JSON aliases linking features between chromosomes/genomes
* Fails if JSON invalid
*/
process generateFeaturesJSON {
  tag{tag}
  label 'json'

  input:
    set val(map), file(pep) from localPepSeqs4Features.mix(remotePepSeqs4Features)
    //WARNINIG! HC OUTPUT CURRENTLY OVERWRITEN BY LC OUTPUT (ORE VICE-VERSA)

  output:
    file "*.json"

  script:
    tag=map.species+"_"+map.version
    // base=map.subgenomes.size() < 2 ? tag : tag+"_"+subgenome
    """
    awk '\$1 ~ /^>/' "${pep}"  | sort -k3,3V | pep2featuresJSON.awk -vname="${tag}" \
    | python -mjson.tool > "${tag}"_annotation.json
    """
}


//REPEAT INPUT FOR EACH SUBGENOME 
localPepSeqs4AliasesRep = Channel.create()
localPepSeqs4Aliases.subscribe onNext: {
  // println it[0]
  if(it[0].containsKey("subgenomes")) {
    for(subgenome in it[0].subgenomes) {
      clone = it[0].clone()
      clone.subgenome = subgenome
      localPepSeqs4AliasesRep << [clone,it[1]]
    }
  } else {
    localPepSeqs4AliasesRep << it
  }
}, onComplete: { localPepSeqs4Aliases.close(); localPepSeqs4AliasesRep.close() }

// process repeatPolyPeps{
  
//   input:  
//     set val(map), file(pep) from localPepSeqs4Aliases
//   output:
//     set val(clone), file(outpep) into localPepSeqs4AliasesRep

//   exec:
//   if(map.containsKey("subgenomes")) {
//       for(subgenome in map.subgenomes) {
//         clone = map.clone()
//         clone.subgenome = subgenome
//         outpep = pep
//         // localPepSeqs4AliasesRep << [clone,it[1]]
//       }
//     // } else {
//     //   localPepSeqs4AliasesRep << it
//     // }
//   }
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


process splitPepSeqsPerSubGenome {
  tag{tag}
  echo true

  input:
    set val(map), file(pep) from localPepSeqs4AliasesRep
    
  output:
    set val(map), file("${tag}.pep") into localPepSeqs4AliasesRepSplit1, localPepSeqs4AliasesRepSplit2
  //   // set val(map), file(pep) optional true into localPepSeqs4AliasesNoSubgenomes

  script:    
    tag=map.species+"_"+map.version+(map.containsKey("subgenome") ? "_"+map.subgenome : "")
    if(map.containsKey("subgenome") && map.subgenomes.size() > 1)  {
       """
       awk '{
        if(\$1 ~ /^>/) {
          split(\$3,arr,":");
          if(arr[3] ~ /${map.subgenome}\$/) {
            current=1; print
          } else {
            current=0
          }
        } else if(current) {
          print
        }
      ;}' pep > ${tag}.pep
      """
    } else {
      """
      cp --preserve=links pep ${tag}.pep
      """
    }
  
  //   //>TRIDC1AG000070.2 pep chromosome:WEWseq_PGSB_20160501:chr1A:409955:411146 gene:TRIDC1AG000070

}


getUniqId = { it.species+it.version+(it.containsKey("subgenome") ? "_"+it.subgenome : "")  }
//COMBINE AND FILTER DUPLICATED CHANNEL TO ALLOW ALL VS ALL DATASETS COMPARISONS
remotePepSeqs4AliasesCombined = remotePepSeqs4Aliases1.mix(localPepSeqs4AliasesRepSplit1).combine(remotePepSeqs4Aliases2.mix(localPepSeqs4AliasesRepSplit2))
  .filter { getUniqId(it[0]) < getUniqId(it[2])  }  //[species,version,file.pep]

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
  tag{tag}
  label 'MMseqs2'
  errorStrategy 'ignore'

  input:
    set val(mapA), file(pepA), val(mapB), file(pepB) from remotePepSeqs4AliasesCombined

  output:
     set val(tagA), val(tagB), file("*.tsv") into pairedProteins

  // when:
  //   !pepA.isEmpty() && !pepB.isEmpty()

  script:
    tagA=mapA.species+"_"+mapA.version //no subgenome spec to be pass to next process as used directly in JSON
    tagB=mapB.species+"_"+mapB.version     
    tag=tagA+(mapA.containsKey("subgenome") ? "_"+mapA.subgenome : "")+"_VS_"+tagB+(mapB.containsKey("subgenome") ? "_"+mapB.subgenome : "")
    // labtag=mapA.toString()+" VS "+mapB.toString()
    """
    mmseqs easy-search ${pepA} ${pepB} ${tag}.tsv \${TMPDIR:-/tmp}/${tag} \
    --greedy-best-hits --threads ${task.cpus} -v 1
    """
}

// //  remotePepSeqs4Aliases1.combine(remotePepSeqs4Aliases2).filter { it[0] != it [3]  && it[0]+it[1] < it[3]+it[4]} .subscribe { println  "NOPE: $it" }

/*
* Generate JSON aliases linking features between chromosomes/genomes
*/
process generateAliasesJSON {
  tag{tag}
  label 'json'

  input:
    set(val(tag1), val(tag2), file(paired)) from pairedProteins

  output:
    file "*.json"

  script:
    tag=tag1+"_VS_"+tag2
    namespace1=tag1+":"+tag1+"_annotation"
    namespace2=tag2+":"+tag2+"_annotation"
    """
    blasttab2json.awk -vnamespace1=${namespace1} -vnamespace2=${namespace2} ${paired} \
    | python -mjson.tool > ${tag}_aliases.json
    """
}

