# Steps required to format some of the included data sets before they can be used as input for the pipeline

## Aegilops tauschii

Mismatch between chr identifiers in reference fasta and gff. Make them consistent 1D,2D,...,7D

```
awk '/^Chr/' annotation/AET_High_confidence_gene.gff3  | sed -E 's/^(Chr)([0-7])/\2D/' > annotation/AET_High_confidence_gene_mod.gff3
awk '/^Chr/' annotation/AET_Low_confidence_gene.gff3  | sed -E 's/^(Chr)([0-7])/\2D/' > annotation/AET_Low_confidence_gene_mod.gff3
sed -E 's/^>(.*chromosome )([1-7]D)(.*)$/>\2 \1\2\3/' AtGSP.fa > AtGSP_mod.fa
```

## Campala *Lr22a*

This is a case where all required information was available in one place, and the approach for extracting that information may be generalised if additional, similarly-formated data sets were to be included.

Download assembly with annotation, and parse:

```
wget -O LS480641.1.embl "http://www.ebi.ac.uk/ena/data/view/LS480641.1&display=text"
python bin/emblparse_Campala_Lr22a.py --infile LS480641.1.embl --outfile LS480641.1.aa.fasta
```

Create dummy "index" file `LS480641.1.len` with a single line:

```
2D  563502314
```

Get chromosome FASTA

python bin/embl_2_fasta.py < LS480641.1.embl | sed '1 s/^.*$/>2D/' >  LS480641.1.fasta


### Svevo *Triticum turgidum ssp. durum*

```bash
for s in $(seq 934111 934124); do
  time wget -O LT${s}.1.embl "https://www.ebi.ac.uk/ena/data/view/LT${s}&display=text";
done
#real    436m21.530s
```

```bash
time paste \
  <(ls LT9341??.1.embl) \
  <( for c in $(seq 1 7); do for g in A B; do echo $c$g; done; done) \
  | while read f c;
    do
    ~/.nextflow/assets/plantinformatics/pretzel-input-generator/bin/emblparse.py \
      -i $f -o ${f%.embl}.aa.fasta -a Svevo -c $c;
    done
```

Collect all protein sequences in one file and record pseudochromosome lengths in another.

```
cat LT9341??.1.aa.fasta  > Svevo.aa.fa
grep -hm2 -e '^RP' -e 'chromosome'  *.embl  | grep -oE -e '[0-9][AB]' -e '[0-9]{2,}' | paste - -  > Svevo.len
```

Extract assembly FASTA

```
for f in Svevo/ENA/*.embl; do ~/.nextflow/assets/plantinformatics/pretzel-input-generator/bin/embl_2_fast
a.py < ${f}; done | sed -E 's/^>.*(..)$/>\1/' > Svevo/ENA/Svevo.fasta
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

Fix chromosome names in FASTA

sed -i 's/^>Tu\([1-7]\)/>\1A/' WheatTu.genome.fasta

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
paste - - < TRIDC_WEWseq_PGSB_20160501_Proteins_HighConf_REPR_on_WEW2.0.fasta | grep -v -e scaffold -e ':\*' | tr '\t' '\n' > TRIDC_WEWseq_PGSB_20160501_Proteins_HighConf_REPR_on_WEW2.0_chromosomes.fasta
```

## Thinopyrum elongatum

```sh
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/799/875/GCA_011799875.1_ASM1179987v1/GCA_011799875.1_ASM11
79987v1_genomic.fna.gz
time pigz -dcp2 GCA_011799875.1_ASM1179987v1_genomic.fna.gz \
| fasta_formatter \
| sed -re 's/^>.*(scaffold[0-9]+.*|[1-7]E),.*/>\1/' \
> Thinopyrum_elongatum_D-3458_chromosomes.fa
```

Try again - annotated data

```sh
time wget ftp://download.big.ac.cn/gwh/Plants/Thinopyrum_elongatum_Thinopyrum_elongatum-REFERENCE-SDAU-1.0.fa_GWHABKY00000000/GWHABKY00000000.genome.fasta.gz ftp://download.big.ac.cn/gwh/Plants/Thinopyrum_elongatum_Thinopyrum_elongatum-REFERENCE-SDAU-1.0.fa_GWHABKY00000000/GWHABKY00000000.gff.gz ftp://download.big.ac.cn/gwh/Plants/Thinopyrum_elongatum_Thinopyrum_elongatum-REFERENCE-SDAU-1.0.fa_GWHABKY00000000/GWHABKY00000000.feature.gz ftp://download.big.ac.cn/gwh/Plants/Thinopyrum_elongatum_Thinopyrum_elongatum-REFERENCE-SDAU-1.0.fa_GWHABKY00000000/GWHABKY00000000.RNA.fasta.gz ftp://download.big.ac.cn/gwh/Plants/Thinopyrum_elongatum_Thinopyrum_elongatum-REFERENCE-SDAU-1.0.fa_GWHABKY00000000/GWHABKY00000000.CDS.fasta.gz ftp://download.big.ac.cn/gwh/Plants/Thinopyrum_elongatum_Thinopyrum_elongatum-REFERENCE-SDAU-1.0.fa_GWHABKY00000000/GWHABKY00000000.Protein.faa.gz

#Update IDs
sed -Ei 's/GWHABKY0000000([1-7])/\1E/' GWHABKY00000000.Protein.faa
sed -Ei 's/>GWHABKY0000000([1-7])/>\1E/' GWHABKY00000000.genome.fasta
sed -Ei 's/GWHABKY0000000([1-7])/\1E/' GWHABKY00000000.gff #may not need it

#Re-format peps to Ensembl-like
bin/bigdPep_2_ensemblPep.awk GWHABKY00000000.Protein.faa > GWHABKY00000000.Protein.ensembl.faa
```


## Oryza sativa

Until EP datasets use is put in sync with local, this could be a way to ensure IRGSP annotations are processed correctly

```
fasta_formatter < Oryza_sativa.IRGSP-1.0.pep.all.fa | paste - - | filterForRepresentative.awk > Oryza_sativa.IRGSP-1.0.pep.REPR.fa
```



## Triticum aestivum IWGSC v2

Lacking annotations, transplant from v1. Get genic sequences from v1:

```
< iwgsc_refseqv1.0_HighConf_2017Mar13.gff3 | awk -vFS="\t" -vOFS="\t" '$3=="gene"{split($9,a,";");split(a[1],b,"=");print $1":"$4"-"$5; print b[2] > "transplant/HC.ids"}'  | xargs samtools faidx --length 10000000 161010_Chinese_Spring_v1.0_pseudomolecules.fasta  > transplant/retrieved_HC.fa
```

## Legacy datasets

### Triticum aestivum IWGSC CSS v2

Source files from https://urgi.versailles.inra.fr/download/iwgsc/Survey_sequence/

```
for F in *fa.gz; do base=${F##*/}; prefix=${base%%-*}; zcat $F | sed -E "s/>([0-9]+).*/>${prefix/_v2/}_\1/"; done > contigs/IWGSC_CSSv2.fa
```

Gene predictions - extract genomic sequences (including with introns) but no additional  up/down stream sequence to preserve coordinates.

```
awk -vOFS="\t" 'NR>1{split($3,arr,"_"); gsub(/\[[-+]\]/,"",$4); print $1,toupper(arr[3])"_"arr[5]":"$4}' survey_sequence_gene_models_MIPS_v2.2_Jul2014/ta_IWGSC_MIPSv2.2_HighConf_REPR_2014Jul18.tab | tee models.info | cut -f2 | xargs samtools faidx ../contigs/IWGSC_CSSv2.fa | fasta_formatter | paste - - > genic.tmp
paste models.info genic.tmp | awk '{gsub(/>/,"",$3); if($2==$3)print ">"$1"\n"$4}' > ta_IWGSC_MIPSv2.2_HighConf_REPR_genic.fasta
```



### Triticum aestivum IWGSC CSS v3


## Ad-hoc processing of pre-computed gffs

```
bin/gff_2_pretzel.py \
  --infile markers/bristol_SNP.summary.gff \
  --name Triticum_aestivum_IWGSC_RefSeq_v1.0_bristolSNP_markers \
  --namespace bristolSNP \
  --short-name bristolSNP \
  --parent Triticum_aestivum_IWGSC_RefSeq_v1.0 \
  --source 'https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.0/' \
  --citation https://doi.org/10.1126/science.aar7191 | pigz -11cp2 > markers/bristol_SNP.json.gz
```

Selected subset

```
for F in {DArT,axiom_820k,iSelect,bristol_SNP}*.summary.gff; do
  NAME=${F%*.summary.gff}
  echo ${NAME}
  time ./gff_2_pretzel.py \
  --infile ${F} \
  --name Triticum_aestivum_IWGSC_RefSeq_v1.0_${NAME}_markers \
  --namespace ${NAME} \
  --short-name ${NAME} \
  --parent Triticum_aestivum_IWGSC_RefSeq_v1.0 \
  --source 'https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.0/' \
  --citation https://doi.org/10.1126/science.aar7191 | pigz -9cp2 > ${NAME}_markers.json.gz
done
```

Axiom 35k is a subset of 820k (https://www.cerealsdb.uk.net/cerealgenomics/CerealsDB/)

```
F=axiom_820k.summary.gff
NAME=axiom_35k
time fgrep -wf <(cut -f8 35k_probe_set_IWGSCv1.tsv | tail -n+2) axiom_820k.summary.gff \
| ./gff_2_pretzel.py \
  --name Triticum_aestivum_IWGSC_RefSeq_v1.0_${NAME}_markers \
  --namespace ${NAME} \
  --short-name ${NAME} \
  --parent Triticum_aestivum_IWGSC_RefSeq_v1.0 \
  --source 'https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.0/' \
  --citation https://doi.org/10.1126/science.aar7191 | pigz -9cp2 > ${NAME}_markers.json.gz
```


## Sugracane & sorghum


```sh
mv \
  local/R570_V1/Saccharum_officinarum_X_spontaneum_var_R570.mainGenome.fasta \
  local/R570_V1/Saccharum_officinarum_X_spontaneum_var_R570.mainGenome.fasta.BAK
sed \
  -e 's/SB_//g' \
  -e 's/chrom/chr/g' \
  -e 's/recombinant/rec/g' \
  -e 's/alternative/alt/g' \
  -e 's/primary/prim/g' \
  local/R570_V1/Saccharum_officinarum_X_spontaneum_var_R570.mainGenome.fasta.BAK \
  | tr '|:' '__' \
  > local/R570_V1/Saccharum_officinarum_X_spontaneum_var_R570.mainGenome.fasta
```

```sh
awk 'NR>1{gsub("/[ACGT]","",$2);print ">"$1"\n"$2}' \
  local/markers/Axiom48k_probes.txt \
  | tr -d '[]' \
  > local/markers/Axiom_48k_probes.fa
```

Not used:

```sh
~/.nextflow/assets/plantinformatics/pretzel-input-generator/bin/gff_2_pretzel.py \
  --infile local/markers/SNPs_vs_sorghum_v2.1_fmtsix.gff \
  --name Sorghum_bicolor_NCBIv3_AX \
  --namespace AX \
  --short-name AX \
  --parent Sorghum_bicolor_NCBIv3 \
  | pigz -11cp2 > results/JSON/Sorghum_bicolor_NCBIv3_AX.json.gz
```


### Sugarcane markers-contigs

It appears that 48k markers are not specific enough - no surprise given probe lengths 
Given the provided contig sequences and a file mapping mmarker IDs to contig IDs,
can we use the contigs as "extended marker sequences?

```sh
head $MAP 
probeset_id     affy_snp_id     cust_id contig  BestProbeset    dbSNP.RS.ID     Chromosome      Physical.Position       Strand  Flank   Allele.A        Allele.B        classification  source
AX-117128103    Affx-115558822  S_126739127     step3_rep_c58742        1       ---     ---     ---     +       TTGTAGTAATAATATTCAAAAAGAAAAGCAAAAGT[C/T]AAGAGGAGAACAGAGGTGCAGCATTCATGGAACTT     T       C       secondchoice    get_variants_for_offtarget.exclude_ontarget.submission
AX-117128112    Affx-115558813  S_5219223       step3_c1924     1       ---     ---     ---     +       GCTGCTTTTCTTCCTCACCGTACCTCTCAAATCGT[C/T]CATGTATTGTTTTTGTAAGGCCCCGTTTTGATCCT     T       C       preferred       class2.submission
AX-117128115    Affx-115558810  S_13564454      step3_c4483     1       ---     ---     ---     +       GCATCCTATCTTGAATCTGTTAGTTGAAGTTCAGG[C/T]TTCACTGTCAGGCCACTTTCATAGTTGTCTTCAGC     T       C       preferred       class2+3.submission
AX-117128120    Affx-115558805  3431655 NA      1       ---     ---     ---     +       CATGCTTACTGAGTCGTGCGAAGACATCACTGGAT[C/G]ACAAAAAAAA      C       G       Infinium        Infinium
AX-117128122    Affx-115558803  S_167692905     step3_rep_c76251        1       ---     ---     ---     +       TATATGTTTTTGAGGGAGCCTCGCCAAACTTTTAG[A/G]GTCTCCTGCAAACAGTGGTTGATGAGAGTAGTCTT     A       G       preferred       class2.submission
AX-117128126    Affx-115558799  S_76015262      step3_rep_c29698        1       ---     ---     ---     +       CTCCAGAATTTTGCTTCAGAAATAGCCTTCGCAAC[C/T]AGCATTCATGCCTGGCTAGCCCCTGCACTGTATGG     T       C       preferred       class1.submission
AX-117128129    Affx-115558796  S_66101515      step3_rep_c25171        1       ---     ---     ---     +       GATGGACCACAATATAGGTAGTAAATGAAACGTGA[A/G]GGAGTGCAGGACTGCAGTGGATCCATC     A       G       secondchoice    get_variants_for_offtarget.exclude_ontarget.submission
AX-117128131    Affx-115558794  S_93811100      step3_rep_c40457        1       ---     ---     ---     +       ACTCGCAGATTAGCGTAACCCCGAAGAGAAAGAGG[A/G]CAAGTCATATGCTCAGAGTATCCAAAAAGAGTAAT     A       G       preferred       class2.submission
AX-117128137    Affx-115558788  S_24992350      step3_c974      1       ---     ---     ---     +       TGGTGAGAGGAAAGATTTGGCATTACAACTGGAGT[C/T]TAAGTCTGCACAAGTTGATGAGCTTCAGCTGATGG     T       C       preferred       class2+3.submission
```

We have one-to-many relatioships between marker IDs and contig IDs


```sh

 cut -f1 $MAP | sort | uniq -c | sort -nr | head
      1 probeset_id
      1 AX-94584309
      1 AX-91428614
      1 AX-91400098
      1 AX-91385174
      1 AX-91192353
      1 AX-91047534
      1 AX-90934140
      1 AX-90637837
      1 AX-90595263

cut -f4 $MAP | sort | uniq -c | sort -nr | head
    866 NA
     20 step3_rep_c19975
     18 step3_c21649
     17 step3_rep_c17089
     17 step3_rep_c12795
     16 step3_rep_c75598
     16 step3_rep_c38807
     16 step3_rep_c25975
     16 step3_rep_c21988
     15 step3_rep_c40546
```

Additionally, we have no contig ID for 866 of 47802 markers.


```sh
awk  'BEGIN{FS=OFS="\t"}; NR==FNR && NR>1{IDS[$4]=$1};NR!=FNR{gsub(">","",$1); print $1,IDS[$1]}' $MAP <(paste - - < $CONTIGS )| head

awk  'BEGIN{FS="\t"; OFS="\n"}; 
  NR==FNR && NR>1 { 
    gsub(">","",$1); 
    ID2FA[$1]=$2
  };
  NR != FNR {
    if($4 in ID2FA) {
      print ">"$1,ID2FA[$4];
    }
  }
' <(paste - - < $CONTIGS ) $MAP > local/markers/Axiom48k_contigs.fa 


```

Currently unused addition to conf/sugarcane.config


```
      [
        name: 'Axiom48k',
        // fasta: 'local/markers/NGS_000105_contigs_new_ids.fa',  //local or remote, either gz or not
        fasta: 'local/markers/Axiom48k_contigs.fa',  //local or remote, either gz or not
        seqtype: 'genomic', //markers|transcripts|cds|genomic
        target: [
          [ species: 'Saccharum_officinarum_X_spontaneum', version: 'R570_V2' ]
        ]
      ]
    ]
```

