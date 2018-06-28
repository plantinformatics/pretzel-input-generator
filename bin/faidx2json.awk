#!/usr/bin/awk -f

BEGIN {
  print "{";
  OFS="    "; #used for indent only
  z="";
  print "  \"name\": \""name"\",";
  print "  \"blocks\": [";
}
{
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
