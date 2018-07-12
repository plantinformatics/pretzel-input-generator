#!/usr/bin/awk -f

BEGIN {
  IGNORECASE=1; #for regex matches
  OFS="    "; #used for indent only
  z="";
  print "{";  
  print "  \"name\": \""name"\",";
  print "  \"blocks\": [";
}
{
  gsub(/chr_?/,"",$1);
  if(NR>1) 
    print z,"},";
  print z,"{";
  print z,z,"\"featureType\": \"linear\",";
  print z,z,"\"range\": [";
  print z,z,z,"1,";
  print z,z,z,$2;
  print z,z,"],";
  print z,z,"\"scope\": \""$1"\"";
}
END {
  print z,"}";
  print "  ]";
  print "}";
}
