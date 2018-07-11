#!/usr/bin/awk -f

BEGIN {
  FS="\t";
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
  if($3 =="mRNA") {
    gsub("\"","");
    split($9,arr,";| "); 
    for(i in arr) {
      split(arr[i], pair, "=");      
      if(pair[1]=="ID") {
        transcript=pair[2];
      } else if(pair[1]=="Parent") {
        gene=pair[2];
      }
    }
    if(transcript in repr) { 
      print ">"transcript" pep chromosome:"version":"$1":"$4":"$5" gene:"gene"\n"repr[transcript]; 
    }
  }
}