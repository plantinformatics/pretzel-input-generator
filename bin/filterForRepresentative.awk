#!/usr/bin/awk -f

BEGIN {
  FS = "\t";
  OFS = "\t";
}

{
  split($1,arr," "); #GET GENE FIELD
  split(arr[4],gene,":"); #GET GENE FIELD
  ID=gene[2]; #GET GENE ID
  sub(/^>[^ ]+/, ">"ID); #USE GENE ID AS FASTA IDENTIFIER (NOT TRANSCRIPT ID)
  if(!(ID in storedIDs) || length($2) > length(StoredSeqLines[ID])) { #FIRST OCCURANCE OT LONGER THAN STORED
    storedIdLInes[ID] = $1;
    # print "storing "$1
    storedSeqLines[ID] = $2;
    # print "storing "$2
  }
}

END {
  for (key in storedIdLInes) {
    print storedIdLInes[key];
    print storedSeqLines[key];
  }
}
