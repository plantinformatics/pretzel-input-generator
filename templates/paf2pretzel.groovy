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
    def arr = line.split('\t')
    //NOTE!!! 0-based coordinates, half-open i.e. first position included, last position excluded, or "<a,b)"
    (QNAME, QLEN, QSTART, QEND, STRAND, TNAME, TLEN, TSTART, TEND, MATCHES, ALNLEN, MAPQ) = arr[0..11] //First 12 are positional fields
    def TAGS = arr[12..-1] //variable number of key-value tag pairs

    // int qlen = QLEN.toInteger()
    int qstart = QSTART.toInteger()
    int qend = QEND.toInteger()
    int tstart = TSTART.toInteger()
    int tend = TEND.toInteger()

    int matches = MATCHES.toInteger()
    // int alnlen = ALNLEN.toInteger()
    // double aligned_identity = matches / ALNLEN.toInteger()
    double query_identity = matches / QLEN.toInteger()
    double query_coverage =  (qend-qstart)/ QLEN.toInteger()
    // int mapq = MAPQ.toInteger()

    if(query_identity >= min_identity) {
      def kosher = true;
      if(query_identity < 1) { //Not a 100% match, check if any MM in last 3 bases
        TAGS.each { tag ->
          if(tag.startsWith('cs:Z')) {
            if(STRAND == '-' && (tag =~ /^cs:Z:=[acgtnACGTN]{3,}/).count == 0){
              // println "minus 3'\t"+line
              kosher = false;
            } else if(STRAND == '+' && (tag.substring(tag.lastIndexOf('=')+1).matches('^[acgtnACGTN]{3,}\$')) == false) {
              // println "plus 3'\t"+line
              kosher = false;
            }
          }
        }
      }
      if(kosher) {
        def key = TNAME.replaceFirst("^(C|c)(H|h)(R|r)[_]?","")
        if(!scope.containsKey(key)) {
          scope << [(key) : []]
        }
        scope[key] << [
          "name" : QNAME,
          "value" : [ tstart+1, tend ],
          "identity" : query_identity,
          "coverage" : query_coverage
        ]
      }
    }
    // scope[key] << ["name" : QNAME, "value" : [ alnStart, alnEnd ]]
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
('gzip ${tag}_markers.json'.execute()).waitFor()
