#!/usr/bin/env python
from __future__ import print_function
import argparse, sys, os
from Bio import SeqIO

def eprint(*args, **kwargs):
  print(*args, file=sys.stderr, **kwargs)

parser = argparse.ArgumentParser(prog=os.path.basename(__file__),
  formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=60),
# parser = argparse.ArgumentParser(os.path.basename(__file__),
  description="Parse ENA embl e.g. from the output of wget -O LS480641.1.embl \"http://www.ebi.ac.uk/ena/data/view/LS480641.1&display=text\"")
parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'),
  default='sys.stdin', help='provide input file name or use stdin')
parser.add_argument('-o', '--outfile', nargs='?', type=argparse.FileType('w'),
  default='sys.stout', help='provide output file name or use stout')
parser.add_argument('-a', '--assembly', help='Assembly version/name', type=str, required=True)
parser.add_argument('-c', '--chromosome', help='Chromosome name', type=str, required=True)


args = parser.parse_args()
f = args.infile
of = args.outfile

# eprint(args.chromosome)
# eprint(args.assembly)
# exit(0)

for record in SeqIO.parse(f, "embl"):
  eprint(record.id)
  eprint("Record {} has {} features".format(record.id, len(record.features)))
  for feature in record.features:
    if feature.type == 'gene':
      geneId = feature.qualifiers['locus_tag'][0]
      coords = feature.location

      # eprint(coords)
      # eprint(feature)

    if feature.type == 'mRNA':
      transcriptId = feature.qualifiers['locus_tag'][0]
      # eprint(feature)
      representative = True if 'is_repr=1' in feature.qualifiers['note'] else False
      # eprint("repr:")
      # eprint(representative)
      # exit(0)
    if feature.type == 'CDS' and representative and 'translation' in feature.qualifiers:
      # eprint(feature)
      aa = feature.qualifiers['translation'][0]
      of.write(">{0} pep chromosome:{1}:{2}:{3}:{4} gene:{5}\n".format(transcriptId,args.assembly,args.chromosome,coords.start,coords.end,geneId))
      of.write("{}\n".format(aa))

f.close()
of.close()
      # break

#TRIDC1AG000120.1 pep chromosome:WEWseq_PGSB_20160501:1A:616966:619040 gene:TRIDC1AG000120

#SeqFeature(CompoundLocation([FeatureLocation(ExactPosition(25975396), ExactPosition(25975637), strand=-1), FeatureLocation(ExactPosition(25970652), ExactPosition(25971464), strand=-1)], 'join'), type='CDS', location_operator='join')