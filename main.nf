#!/usr/bin/env nextflow
//echo true

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
  tag{tag}
  
  input:
    val species from assemblySpeciesList
    val version from assemblyVersionList
    val urlprefix
    val pepsuffix
    val idxsuffix

  output:
    set val(species), val(version), file(idx) into remoteIndices
    set val(species), val(version), file(pep) into remotePepSeqs //remotePepSeqs4Features, remotePepSeqs4Aliases1, remotePepSeqs4Aliases2
    // file 'pep*' into remotePepSeqsNoLabs
 
  script:
    tag=species+"_"+version
    idxurl=urlprefix+eprelease+"/fasta/"+species.toLowerCase()+"/dna_index/"+species+"."+version+idxsuffix
    pepurl=urlprefix+eprelease+"/fasta/"+species.toLowerCase()+"/pep/"+species+"."+version+pepsuffix
    """
    curl $idxurl > idx
    curl $pepurl | gunzip --stdout | head -20000 > pep #"${tag}.pep"
    """
}

/*
* Only keep "representative" splice form for each gene, typically suffix ".1", some times "-01"
*/
process filterForRepresentativePeps {
  tag{tag}
  //label 'fastx'
  input:
    set val(species), val(version), file(pep) from remotePepSeqs
  
  output:
    set val(species), val(version), file("${tag}.pep") into remotePepSeqs4Features, remotePepSeqs4Aliases1, remotePepSeqs4Aliases2

  script:
  tag=species+"_"+version
    """
    cat ${pep} | filterForRepresentative.awk > "${tag}.pep"
    """
}

/*
* Generate genome blocks definitions JSON for pretzel 
* Fails if JSON invalid
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
    | faidx2json.awk -vname="${tag}" \
    | python -mjson.tool > "${tag}"_genome.json
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
    set val(species), val(version), file(pep) from remotePepSeqs4Features
  
  output:
    file "*.json"

  script:    
    tag=species+"_"+version
    """
    awk '\$1 ~ /^>/' "${pep}"  | sort -k3,3V | pep2featuresJSON.awk -vname="${tag}" \
    | python -mjson.tool > "${tag}"_annotation.json
    """
}




/*
* Identify best hit for each pep
*/
process pairProteins {
  tag{tag}
  label 'MMseqs2'

  input:
    each one from remotePepSeqs4Aliases1.toList()
    each another from remotePepSeqs4Aliases2.toList()
  
  output:
    set val(tag1), val(tag2), file("*.tsv") into pairedProteins

  when:
    one != another && one[0]+one[1] < another[0]+another[1] //DON'T COMPARE WITH ITSELF, DO ONE-WAY ONLY BASED ON LEX ORDER

  script:
    tag1=one[0]+"_"+one[1]
    tag2=another[0]+"_"+another[1]
    tag=tag1+"_VS_"+tag2
    """
     # for f in ${one[2]} ${another[2]}; do
    #   awk '{if(\$1 ~ /^>/){split(\$4,gene,":"); print ">"gene[2]}else{print}' \${f} > \${f##*/}
    # done
    mmseqs easy-search ${one[2]} ${another[2]} ${tag}.tsv \${TMPDIR:-/tmp} --greedy-best-hits
    """
}

/*
* Generate JSON aliases linking features between chromosomes/genomes
*/
process generateAliasesJSON {
  tag{tag}
  label 'json'

  input:
    set val(tag1), val(tag2), file(paired) from pairedProteins

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


