#!/usr/bin/env groovy

import static groovy.json.JsonOutput.*

buscos = new File('${tsv}').text
out = new File('${tag}.json')

//RECORD chromosome/scaffold sizes
idx = new File('${fai}').text
lengths = [:]
idx.eachLine { line ->
  if(line.toLowerCase() =~ /^(chr|[0-9]|x|y|hic_scaffold)/ ) {
toks = line.split('\t')
lengths.(toks[0].replaceFirst("^(C|c)(H|h)(R|r)[_]?","")) = toks[1].toInteger()
  }
}
//AGGREGATE DATA MAP
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
//current = [ "scope": k, "featureType": "linear", "features": []]
current = [ "scope": k, "featureType": "linear", "range": [1, lengths[k]], "features": []]
features.each { feature ->
  current.features << feature
}
annotation.blocks << current
  }
}
out.text = prettyPrint(toJson(annotation))
'gzip ${tag}.json'.execute()