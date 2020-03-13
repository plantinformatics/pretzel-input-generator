#!/usr/bin/gawk -f

#example usage: NS="SSR"; ./tsvPairs2aliasesJson.awk -vnamespace1=${NS} -vnamespace2=${NS}  two_column_pairs.tsv > aliases.json

BEGIN {
  print "[";
  OFS="    "; #used for indent only
  z="";
}
{
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
