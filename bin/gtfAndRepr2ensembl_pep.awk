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
  if($3=="transcript") {
    gsub("\"","");
    split($9,arr,";| "); 
    for(i in arr) {
      if(arr[i]=="gene_id") {
        gene=arr[i+1];
      } else if (arr[i]=="transcript_id") {
        transcript=arr[i+1];
      }
    }
    if(transcript in repr) { 
      print ">"transcript" pep chromosome:"version":"$1":"$4":"$5" gene:"gene"\n"repr[transcript]; 
    }
  }
}
