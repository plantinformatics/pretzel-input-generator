#!/usr/bin/env groovy

import static groovy.json.JsonOutput.*

pafContent = new File('${paf}').text
out = new File('${tag}_markers.json')
counts = new File('${tag}_markers.counts')
min_identity = 0.9

//AGGREGATE DATA MAP
def annotation = [:]
annotation.meta = [:]

annotation.'public' = true
annotation.name = "${tag}_markers"
annotation.namespace = "${tag}"
annotation.parent = "${genome}"
annotation.meta.shortName = "${meta.markers.name}"

TreeMap scope = [:] //keep keys sorted as the corresponding blocks get displayed in order in pretzel
pafContent.eachLine { line ->
  if(!line.startsWith('@')) {

    def arr = line.split('\t')
    //NOTE!!! 0-based coordinates, half-open i.e. first position included, last position excluded, or "<a,b)"
    (QNAME, QLEN, QSTART, QEND, STRAND, TNAME, TLEN, TSTART, TEND, MATCHES, ALNLEN, MAPQ) = arr[0..11] //First 12 are positional fields
    def TAGS = arr[12..-1] //variable number of key-value tag pairs

    // int qlen = QLEN.toInteger()
    // int qstart = QSTART.toInteger() + 1
    // int qend = QEND.toInteger()
    int tstart = TSTART.toInteger() + 1
    int tend = TEND.toInteger()

    // int matches = MATCHES.toInteger()
    // int alnlen = ALNLEN.toInteger()
    double identity = MATCHES.toInteger() / ALNLEN.toInteger()
    // int mapq = MAPQ.toInteger()

    if(identity >= min_identity) {
      def key = TNAME.replaceFirst("^(C|c)(H|h)(R|r)[_]?","")
      if(!scope.containsKey(key)) {
        scope << [(key) : []]
      }
      scope[key] << [
        "name" : QNAME,
        "value" : [ tstart, tend ],
        "identity" : matches/alnlen
      ]
    }
    // scope[key] << ["name" : QNAME, "value" : [ alnStart, alnEnd ]]
  }
}
//GROUP TOGETHER FEATURES FROM/IN SAME BLOCK
annotation.blocks = []
scope.each { k, features ->
  current = [ "scope": k, "featureType": "linear", "features": []]
  // // current = [ "scope": k, "featureType": "linear", "range": [1, lengths[k]], "features": []]
  current.features = features
  annotation.blocks << current
}
//RECORD NUM FEATURS PER BLOCK
counts.withWriterAppend{ wr ->
  annotation.blocks.each {
    wr.println annotation.name+"\\t"+it.scope+"\\t"+it.features.size()
  }
}
out.text = prettyPrint(toJson(annotation))
'gzip ${tag}_markers.json'.execute()
