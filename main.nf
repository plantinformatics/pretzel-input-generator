#!/usr/bin/env nextflow
// echo true

//INPUT PARAMS
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

// localIds1.merge(Channel.fromPath(params.localGtf)).merge(Channel.fromPath(params.localPep)).set{ localInputGtfPep }
// localIds2.merge(Channel.fromPath(params.localIdx)).set { localIndices } 

//ARRANGE INPUTS FOR PROCESSES
localInputGtfPep = Channel.create()
localIndices = Channel.create()
params.localAssembly.each {
  //EXPECT TO HAVE SOME DATASETS WITH gff3 instead
  if(it.containsKey("gtf")) {
    localInputGtfPep << [it,file(it.gtf),file(it.pep)]
  }
  //ALL SHOULD HAVE AN INDEX
  if(it.containsKey("idx")) { 
    localIndices << [it,file(it.idx)]
  }
}
localInputGtfPep.close()
localIndices.close()

// localInputGtfPep.subscribe {
//    println ("HERE: $it")
// }

/*
* Generate genome blocks definitions JSON for pretzel 
* Fails if JSON invalid
*/
process generateGenomeBlocksJSON {
  tag{tag}
  label 'json'

  input:
    set val(map), file(idx) from localIndices//.mix(remoteIndices)

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
* Given a FASTA with representative peps and the corresponding GTF
* output FASTA with representative peps and definition lines 
* mimicking the ensembl plants (EP) format for such data - this can then 
* be piped into the same processes which we use for chewing through EP data
*/
process convertReprFasta2EnsemblPep {
  tag{tag}
  label 'fastx'

  input:
    set (val(map), file(gtf), file(reprPep)) from localInputGtfPep

  output:
    set val(map), file(pep) into localPepSeqs4Features, localPepSeqs4Aliases  
  script:    
    tag=map.species+"_"+map.version
    """
    fasta_formatter < ${reprPep} | gtfAndRepr2ensembl_pep.awk -vversion="${map.version}" -vFS="\\t" - ${gtf} > pep
    """
}



/*
* Download peptide seqs and assembly index files from Ensembl plants
*/
process fetchRemoteDataFromEnsemblPlants {
  tag{map}
  
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
    """
    curl $idxurl > idx
    curl $pepurl | gunzip --stdout | head -20000 > pep
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
    // base=map.subGenomes.size() < 2 ? tag : tag+"_"+subgenome
    """
    awk '\$1 ~ /^>/' "${pep}"  | sort -k3,3V | pep2featuresJSON.awk -vname="${tag}" \
    | python -mjson.tool > "${tag}"_annotation.json
    """
}

// // subgenomes = params.localSubGenomes.tokenize(",")
// // process splitPepSeqsPerSubGenome {
// //   tag{tag}

// //   input:
// //     set val(species), val(version), file(pep) from localPepSeqs4Aliases
// //     each subgenome from subgenomes
// //   output:
// //     set val(species), val(version), file(subpep) into localPepSeqs4AliasesSubgenomes
    
// //   script:
// //   if(subgenomes.size() > 1)  {
// //     tag=species+"_"+version+"_subgenome:"+subgenome
// //     version+="_"+subgenome
// //     """
// //     awk '{if(\$1 ~ /^>/) {split(\$3,arr,":");if(arr[3] ~ /${subgenome}\$/){current=1; print}else{current=0}}else if(current){print};}' pep > subpep
// //     """
// //   } else {
// //     tag=species+"_"+version
// //     """
// //     cp --no-dereference pep subpep
// //     """
// //   }
// //     //>TRIDC1AG000070.2 pep chromosome:WEWseq_PGSB_20160501:chr1A:409955:411146 gene:TRIDC1AG000070

// // }




//COMBINE AND FILTER DUPLICATED CHANNEL TO ALLOW ALL VS ALL DATASETS COMPARISONS
remotePepSeqs4AliasesCombined = remotePepSeqs4Aliases1.combine(remotePepSeqs4Aliases2)
  .filter { it[0].species+it[0].version < it[2].species+it[2].version  }  //[species,version,file.pep]
// remotePepSeqs4AliasesCombined.println()
// [
  // [species:Arabidopsis_thaliana, version:TAIR10], 
  // /home/rad/repos/pretzel-input-generator/work/59/e1b177dcc513678af115b00116e9b6/Arabidopsis_thaliana_TAIR10.pep, 
  // [species:Brachypodium_distachyon, version:v1.0], 
  // /home/rad/repos/pretzel-input-generator/work/82/761d637387c59c8b6ceeb50ea870e2/Brachypodium_distachyon_v1.0.pep]

/*
* Identify best hit for each pep
*/
process pairProteins {
  tag{labtag}
  label 'MMseqs2'

  input:
    set val(mapA), file(pepA), val(mapB), file(pepB) from remotePepSeqs4AliasesCombined 

  output:
     set val(tagA), val(tagB), file("*.tsv") into pairedProteins
  
  script:
    tagA=mapA.species+"_"+mapA.version 
    tagB=mapB.species+"_"+mapB.version 
    tag=tagA+"_VS_"+tagB
    labtag=mapA.toString()+" VS "+mapB.toString()
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

