#!/usr/bin/env groovy

import static groovy.json.JsonOutput.*

samContent = new File('${sam}').text
out = new File('${tag}.json')

//AGGREGATE DATA MAP
def annotation = [:]
annotation.meta = [:]

annotation.name = "${tag}"

TreeMap scope = [:] //keep keys sorted as the corresponding blocks get displayed in order in pretzel
samContent.eachLine { line ->
  if(!line.startsWith('@')) {
    (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL) = line.split('\t')
    int mapq = MAPQ.toInteger()
    int flag = FLAG.toInteger()

    //MULTI-MAPPER?
    boolean isSecondary = ((flag & 256) != 0)
    boolean isSupplementary = ((flag & 2048) != 0)
    boolean isPrimary = (!isSecondary && !isSupplementary)
    // statsKey = isPrimary ? 'primary' : isSecondary ? 'secondary' : isSupplementary ? 'supplementary' : "ERROR - record should be primary or secondary or supplementary"

    if(RNAME != '*') { //Aligned
      // //ORIENTATION
      // boolean reverse = ((flag & 16) != 0)

      //Alignment details
      def cigar = CIGAR
      int alnStart = POS.toInteger() //- getStartClip(cigar)
      int alnEnd = alnStart + getAlignmentLength(cigar)-1 //Get aln len from CIGAR string

      if(!scope.containsKey(RNAME)) {
        scope << [(RNAME) : []]
      }
      scope[RNAME] << ["name" : QNAME, "value" : [ alnStart, alnEnd ]]
    }
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
out.text = prettyPrint(toJson(annotation))


def int getAlignmentLength(String cigar) {
  int len = 0
  StringBuilder segment = new StringBuilder();
  for(int i=0; i < cigar.length(); i++){
    char c=cigar.charAt(i)
    if(Character.isDigit(c)){
      segment.append(c)
    } else {
      // if(c == 'M' || c == '=' || c == 'X' || c == 'D' || c == 'N' || c == 'S' || c =='H') {
      if(c == 'M' || c == '=' || c == 'X' || c == 'D' || c == 'N') {
        len+=segment.toInteger();
      // }else if(c=='I'){
      } else if(c=='I' || c == 'S' || c =='H'){
        //ignore insertions in reference?
        //also ignore soft- and hard-clipping?
      } else if(c=='*'){
        //CIGAR UNAVAILABLE
        return 1 //If no CIGAR in pseudobam then record we may have alnStart but no alnEnd
      } else{
        System.err.println("Unexpected character "+c+" in CIGAR string: "+cigar)
        System.exit(1)
      }
      segment.setLength(0) //clear buffer
    }
  }
  return len;
}

def int getStartClip(String cigar) {
  def m = (cigar =~ /^[0-9]+[HS]/)
  return m.count == 0 ? 0 : (m[0][0..-2]).toInteger()
}

def int getEndClip(String cigar) {
  def m = (cigar =~ /[0-9]+[HS]\$/)
  return m.count == 0 ? 0 : (m[0][0..-2]).toInteger()
}