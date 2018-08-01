#!/usr/bin/awk -f

#example usage: NS="IWGSC_RefSeq_v1.0:IWGSC_RefSeq_v1.0_annotation"; ./blasttab2json.awk -vnamespace1=${NS} -vnamespace2=${NS}  AB_greedybest_HCLC.tsv > AB_greedybest_HCLC.json

BEGIN {
  print "[";
  OFS="    "; #used for indent only
  z="";
}
{
  gsub("\\.[0-9]+$","",$1);
  gsub("\\.[0-9]+$","",$2);
  if(NR>1)
    print z,"},";
  print z,"{";
  print z,z,"\"namespace1\": \""namespace1"\",";
  print z,z,"\"namespace2\": \""namespace2"\",";
  print z,z,"\"string1\": \""$1"\",";
  print z,z,"\"string2\": \""$2"\",";
  #EVIDENCE (OPTIONAL)
  print z,z,"\"evidence\": {";
  print z,z,z,"\"type\": \"alignment\",";
  print z,z,z,"\"tool\": \"mmseqs2\",";
  print z,z,z,"\"identity\": " $3",";
  # print z,z,z,"\"length\": " $4",";
  print z,z,z,"\"coverage1\": " ($8-$7+1)/$13",";  #off by one if pep seq terminated with asterisk !!!
  print z,z,z,"\"coverage2\": " ($10-$9+1)/$14; #",";   
  # print z,z,z,"\"evalue\": " $11","; 
  # print z,z,z,"\"bitscore\": " $12; 
  print z,z,"}";

}
END {
  print z,"}";
  print "]";
}
