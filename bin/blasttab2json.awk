#!/usr/bin/awk -f

#example usage: NS="IWGSC_RefSeq_v1.0:IWGSC_RefSeq_v1.0_annotation"; ./blasttab2json.awk -vnamespace1=${NS} -vnamespace2=${NS} AB_greedybest_HCLC.tsv > AB_greedybest_HCLC.json

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
  print z,z,"\"string2\": \""$2"\"";
}
END {
  print z,"}";
  print "]";
}
