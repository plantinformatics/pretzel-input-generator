# Steps required to format some of the included data sets before they can be used as input for the pipeline

## Aegilops tauschii

Mismatch between chr identifiers in reference fasta and gff. Make them consistent 1D,2D,...,7D

```
awk '/^Chr/' annotation/AET_High_confidence_gene.gff3  | sed -E 's/^(Chr)([0-7])/\2D/' > annotation/AET_High_confidence_gene_mod.gff3
awk '/^Chr/' annotation/AET_Low_confidence_gene.gff3  | sed -E 's/^(Chr)([0-7])/\2D/' > annotation/AET_Low_confidence_gene_mod.gff3
sed -E 's/^>(.*chromosome )([1-7]D)(.*)$/>\2 \1\2\3/' AtGSP.fa > AtGSP_mod.fa
```

## Campala *Lr22a*

This is a case where all required information was available in one place, and the approach for extracting that information may be generalised if additional, similarily-formated data sets were to be included.

Download assembly with annotation, and parse:

```
wget -O LS480641.1.embl "http://www.ebi.ac.uk/ena/data/view/LS480641.1&display=text"
python bin/emblparse_Campala_Lr22a.py --infile LS480641.1.embl --outfile LS480641.1.aa.fasta
```


## *Triticum urartu*

Download assembly and annotation from MBKBase

```
wget http://www.mbkbase.org/Tu/WheatTu.genome.fasta.gz
wget http://www.mbkbase.org/Tu/WheatTu.annotation.tar.gz
tar xzvf WheatTu.annotation.tar.gz
```

Separate high and low confidence genes

```
awk -vOFS="\n" '
NR==FNR {
  if($3 ~ /High/) {
    high[$1]=$2;
  } else {
    low[$1]=$2;
  }
};
NR!=FNR {
  gsub(">","");
  gsub(/\.P[0-9]+/,"",$1);
  if($1 in high) {
    print ">"$1,$2 > "WheatTu.pros.long.HC.fasta"
  } else if($1 in low) {
    print ">"$1,$2 > "WheatTu.pros.long.LC.fasta"
  } else {
    print "Error, expected each gene to be either High or Low confidence!\n"$0; exit 1;
  }
}' WheatTu.gene.evidence <(paste - - < WheatTu.pros.long.fasta)
```

Fix chromosome names in gff, exclude unplaced scaffolds/contigs, remove `.T??` suffixes from gene/transcript name

```
grep -E '^Tu[1-7]' WheatTu.gene.gff | sed -E -e 's/^(Tu)([1-7])/\2A/' -e 's/\.T[0-9]+//g' > WheatTu.gene_mod.gff
```

## Triticum dicoccoides (Wild Emmer) WEW_v2.0


Get data.

```
wget ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/ls/LSYQ01.fasta.gz
```

Decompress and sanitize ids.

```
pigz -dcp2 LSYQ01.fasta.gz | sed -re 's/^>ENA.*(scaffold[0-9]+.*|[1-7][A-B]),.*/>\1/' > WEW_2.0.fasta
```

Align
```
./minimap2-2.16_x64-linux/minimap2 -a -x splice -I 20G -t 20 WEW_2.0.fasta ../Wheat_Wild_Emmer_Zavitan/TRIDC_WEWseq_PGSB_20160501_CDS_HighConf_REPR.fasta > TRIDC_WEWseq_PGSB_20160501_CDS_HighConf_REPR_vs_WEW_2.0.sam
```

* `-I 20G` ensures unsplit index needed for valid header without re-creating the header, an alternative would be to redirect to  `| samtools view -b -T ref.fa - > bam`


```
samtools view -F 2304 TRIDC_WEWseq_PGSB_20160501_CDS_HighConf_REPR_vs_WEW_2.0.sam \
| awk '
  NR==FNR {
    len=0;
    split($6,a,/[[:upper:]]/);
    for(i in a){len+=a[i]};
    loc[">"$1]=$3":"$4":"$4+len
  }
  NR!=FNR {
    if($1 ~ /^>/) {
      ID=$1; ID2=$1;
      sub(/^>/,"",ID2);
      sub(/\.[0-9]+$/,"",ID2);
      print ID,"pep","chromosome:WEW_v2.0:"loc[ID],"gene:"ID2
    } else {
      print
    }
  }
' - <(fasta_formatter < ../Wheat_Wild_Emmer_Zavitan/TRIDC_WEWseq_PGSB_20160501_Proteins_HighConf_REPR.fasta) \
  > TRIDC_WEWseq_PGSB_20160501_Proteins_HighConf_REPR_on_WEW2.0.fasta
```

Now skip genes which are placed on scaffolds not chromosomes

```
paste - - < TRIDC_WEWseq_PGSB_20160501_Proteins_HighConf_REPR_on_WEW2.0.fasta | grep -v -e scaffold -e '*' | tr '\t' '\n' > TRIDC_WEWseq_PGSB_20160501_Proteins_HighConf_REPR_on_WEW2.0_chromosomes.fasta
```