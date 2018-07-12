#!/usr/bin/awk -f

BEGIN {
  IGNORECASE=1; #for regex matches
  print "{";
  OFS="    "; #used for indent only
  z="";
  print "  \"name\": \""name"_genes\","; 
  print "  \"parent\": \""name"\","; 
  print "  \"namespace\": \""name":"name"_annotation\",";
  print "  \"blocks\": [";
}
{
  if($1 ~ /^>/ && $3 ~/^chromosome/) { #>AT3G05780.1 pep chromosome:TAIR10:3:1714941:1719608:-1 gene:AT3G05780    
    split($3,location,":");
    gsub(/chr_?/,"",location[3]);
    split($4,gene,":");
    if(location[3] == currentChromosome) { #ANOTHER GENE IN BLOCK/CHR
      print z,z,z,"},";
    } else {
      if(currentChromosome != "") { #NOT THE FIRST ONE
        print z,z,z,"}";
        print z,z,"],";
        print z,z,"\"featureType\": \"linear\"";
        print z,"},";
      }
      currentChromosome = location[3];
      print z,"{";
      print z,z,"\"scope\": \""currentChromosome"\",";
      print z,z,"\"features\": [";
    }

    print z,z,z,"{";
    print z,z,z,z,"\"range\": [";
    print z,z,z,z,z,location[4]",";
    print z,z,z,z,z,location[5];
    print z,z,z,z,"],";
    print z,z,z,z,"\"name\": \""gene[2]"\"";
  }
}
END {
  print z,z,z,"}";
  print z,z,"],";
  print z,z,"\"featureType\": \"linear\""; 
  print z,"}";
  print "  ]";
  print "}";
}
