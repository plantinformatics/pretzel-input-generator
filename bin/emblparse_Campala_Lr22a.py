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
parser.add_argument('-o', '--outfile', nargs='?', type=argparse.FileType('wa'),
  default='sys.stout', help='provide output file name or use stout')


args = parser.parse_args()
of = args.outfile

for record in SeqIO.parse(args.infile, "embl"):
  eprint(record.id)
  eprint("Record {} has {} features".format(record.id, len(record.features)))
  for feature in record.features:
    if feature.type == 'gene':
      geneId = feature.qualifiers['standard_name'][0]
      coords = feature.location
      # eprint(coords)
      # eprint(feature)
    if feature.type == 'mRNA':
      transcriptId = feature.qualifiers['standard_name'][0]
      # eprint(feature)
    if feature.type == 'CDS':
      aa = feature.qualifiers['translation'][0]
      # eprint(feature)
      of.write(">{0} pep chromosome:{1}:{2}:{3}:{4} gene:{5}\n".format(transcriptId,"Campala_Lr22a","2D",coords.start,coords.end,geneId))
      of.write("{}\n".format(aa))
      # break

#TRIDC1AG000120.1 pep chromosome:WEWseq_PGSB_20160501:1A:616966:619040 gene:TRIDC1AG000120

#SeqFeature(CompoundLocation([FeatureLocation(ExactPosition(25975396), ExactPosition(25975637), strand=-1), FeatureLocation(ExactPosition(25970652), ExactPosition(25971464), strand=-1)], 'join'), type='CDS', location_operator='join')