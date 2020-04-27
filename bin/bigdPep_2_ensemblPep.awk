#!/usr/bin/gawk -f

#Developed for ftp://download.big.ac.cn/gwh/Plants/Thinopyrum_elongatum_Thinopyrum_elongatum-REFERENCE-SDAU-1.0.fa_GWHABKY00000000/GWHABKY00000000.Protein.faa.gz

BEGIN {
  FS="\t";
  version="SDAU-1.0";
}
{
  if($1 ~ /^>/ ) {
    gsub(">","");
    for(i=2; i<=NF; i++) {
      split($i, pair, "=");
      if(pair[1]=="OriID") {
        transcript=pair[2];
      } else if(pair[1]=="OriGeneID") {
        gene=pair[2];
      } else if(pair[1]=="Position") {
        gsub(/ /, "", pair[2])
        n=split(pair[2],arr,":|-|,");
        chr=arr[1]      
        from=arr[2]  
        for(j in arr) {
          if(arr[j] ~ /[0-9]+/) {
            to = arr[j]
          }        
        }
        
      }     
    }
    print ">"transcript" pep chromosome:"version":"chr":"from":"to" gene:"gene; 
  } else {
    print
  }
}
