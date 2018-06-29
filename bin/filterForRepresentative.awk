#!/usr/bin/awk -f

{
  if($0 ~ /^>/) { #ID LINE
    if($1 ~ /\.1$|\-01$/) { #LOOKS LIKE REPRESENTATIVE SPLICE FORM 
      split($4,gene,":"); #GET GENE ID 
      $1=">"gene[2];   #REPLACE SPLICEFORM ID WITH GENE ID
      print;
      repr=1;
    } else {
      repr=0;
    }
  } else if(repr){ #SEQUENCE LINE
    print; 
  }
}
