#!/usr/bin/gawk -f

BEGIN {
  FS="\t";
  IGNORECASE=1; #for regex matches
}
NR==FNR {
  if($1 ~ /^>/ ) {
    gsub(">","");
    gsub(" .*$","");
    key=$1;
  } else {
    repr[key]=$0;
  }
}

NR!=FNR {
  if($3 =="CDS") {
    gsub("\"","");
    split($9,arr,";| ");
    for(i in arr) {
      split(arr[i], pair, "=");
      if(pair[1]=="ID") {
        transcript=pair[2];
      } else if(pair[1]=="Parent") {
        parent=pair[2];
        gene=parent
        gsub(/\.[0-9]+$/,"",gene);
      } else if(pair[1] ~ /^protein(_source)?_id$/)  {
        source=pair[2];
      }
    }
    # print "p="parent,"g="gene,"t="transcript,"s="source
    # if(transcript in repr && !(gene in printed)) { 
    if(!(parent in printed)) {
      #IGNORECASE=1;
      gsub(/chr_?/,"",$1);
      #IGNORECASE=0;
      
      if(source in repr) {
        id = source
      } else if(parent in repr) { #dealing with GFF files being inconsistent...
        id = parent
      } else {
        id = ""
      }
    # if(source in repr && !(gene in printed)) {
      # print ">"transcript" pep chromosome:"version":"$1":"$4":"$5" gene:"gene"\n"repr[transcript];
      if(id) {
        print ">"id" pep chromosome:"version":"$1":"$4":"$5" gene:"gene"\n"repr[id];       
        printed[parent]=1;       
      }
    }
  }
}