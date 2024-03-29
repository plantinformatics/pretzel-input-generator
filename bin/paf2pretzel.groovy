#!/usr/bin/env groovy

import static groovy.json.JsonOutput.*
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream


//@Grab('info.picocli:picocli-groovy:4.1.2') //command line interface
groovy.grape.Grape.grab(group:'info.picocli', module:'picocli-groovy', version:'4.5.0')

@Command(header = [
       //Font Name: Calvin S
       $/@|bold,blue  ╔═╗╔═╗╔═╗  ┌┬┐┌─┐  ╔═╗┬─┐┌─┐┌┬┐┌─┐┌─┐┬   |@/$,
       $/@|bold,blue  ╠═╝╠═╣╠╣    │ │ │  ╠═╝├┬┘├┤  │ ┌─┘├┤ │   |@/$,
       $/@|bold,blue  ╩  ╩ ╩╚     ┴ └─┘  ╩  ┴└─└─┘ ┴ └─┘└─┘┴─┘ |@/$
       ],
       description = "Convert PAF alignment of marker sequences to Pretzel features JSON",
       showDefaultValues = true,
       footerHeading = "%nFootnote(s)%n",
       footer = ["[1] ASCII Art thanks to http://patorjk.com/software/taag/"]
)
@picocli.groovy.PicocliScript
import groovy.transform.Field
import java.security.MessageDigest
import static picocli.CommandLine.*

//should be optional (then generate with blocks, otherwise link to existing dataset)
// @Option(names = ["-f", "--faidx"], description = ["FAI index for the reference FASTA file"], required = true)
// @Field private String fai
@Option(names = ["--min-identity"], description = ["Minimum query sequence identity for sequence placement to be repoted"])
@Field private double minIdentity = 0.9

// @Option(names = ["--min-coverage"], description = ["Minimum query sequence coverage for sequence placement to be repoted"])
// @Field private double minCoverage = 0.7

@Option(names = ["--make-private"], description = ["Make output dataset private (is public by default)"])
@Field private boolean makePrivate = false

@Option(names = ["--parent"], description = ["Parent dataset name"])
@Field private String parent

@Option(names = ["--base-name"], description = ["Basename used for dataset name"])
@Field private String basename

@Option(names = ["--namespace"], description = ["String used for dataset namespace"])
@Field private String namespace

@Option(names = ["--short-name"], description = ["Short display name"])
@Field private String shortName

@Option(names = ["--sequence-type"], description = ["Input sequences type (markers|genomic|transcripts|cds|orf)"])
@Field private String seqType = 'markers'

@Option(names = ["--paf"], description = ["Input PAF file name"])
@Field private String paf = '/dev/stdin'

@Option(names = ["--align-tool"], description = ["Tool used to generate input PAF alignments"])
@Field private String alignTool

@Option(names = ["--align-params"], description = ["Params used to generate input PAF alignments"])
@Field private String alignParams

@Option(names = ["--allowed-target-id-pattern"], description = ["Provide target identifier patter if other than common chromosome naming"])
@Field private String allowedTargetIdPattern

@Option(names = ["-O", "--output"], description = ["JSON output file name"])
@Field private String output = '/dev/stdout'

@Option(names = ["-C", "--out-counts"], description = ["Per-block feature counts report"])
@Field private String outputCounts = '/dev/stderr'

@Option(names = ["-t", "--out-tsv"], description = ["Tab-delimited output"])
@Field private String outputTsv

@Option(names = ["-g", "--out-gff"], description = ["GFF output"])
@Field private String outputGff

@Option(names= ["-h", "--help"], usageHelp=true, description="Show this help message and exit.")
@Field private boolean helpRequested



boolean markerMode = seqType == 'markers'

def pafContent = new File(paf).text
// def out = new File('${tag}_markers.json')
def out = new File(output)
// def counts = new File('${tag}_markers.counts')
def counts = new File(outputCounts)


//AGGREGATE DATA MAP
def annotation = [:]
annotation.meta = [:]

annotation.'public' = !makePrivate
annotation.name = "${basename}_${seqType}"
annotation.namespace = namespace
annotation.parent = parent
annotation.meta.shortName = shortName
annotation.meta.type = seqType.capitalize()
if(alignTool != null && alignParams != null) {
  annotation.meta.evidence = [
    type: "alignment",
    tool: alignTool,
    params: alignParams
  ]
}

String csTag = null
TreeMap scope = [:] //keep keys sorted as the corresponding blocks get displayed in order in pretzel
pafContent.eachLine { line ->
    def arr = line.split('\t')
    //NOTE!!! 0-based coordinates, half-open i.e. first position included, last position excluded, or "<a,b)"
    (QNAME, QLEN, QSTART, QEND, STRAND, TNAME, TLEN, TSTART, TEND, MATCHES, ALNLEN, MAPQ) = arr[0..11] //First 12 are positional fields
    def TAGS = arr[12..-1] //variable number of key-value tag pairs

    // int qlen = QLEN.toInteger()
    int qstart = QSTART.toInteger()
    int qend = QEND.toInteger()
    int tstart = TSTART.toInteger()
    int tend = TEND.toInteger()

    int matches = MATCHES.toInteger()
    // int alnlen = ALNLEN.toInteger()
    // double aligned_identity = matches / ALNLEN.toInteger()
    double query_identity = matches / QLEN.toInteger()
    double query_coverage =  (qend-qstart)/ QLEN.toInteger()
    // int mapq = MAPQ.toInteger()

    // println "${query_identity} >= ${minIdentity} ?"
    // println "${query_identity}\t${query_coverage}\t${query_identity*query_coverage}"
    csTag = null;
    TAGS.each { tag ->
      if(tag.startsWith('cs:Z')) {
        csTag = tag.substring(5)
      }
    }

    // if(query_identity >= minIdentity && query_coverage >= minCoverage) {
    if(query_identity >= minIdentity) {
      def kosher = true;
      //check if TNAME kosher
      if(!((TNAME.toLowerCase() =~ /^(ch|[0-9]{1,2}|x|y|i|v)/) || (TNAME =~ allowedTargetIdPattern) )) {
        kosher = false //don't report placement on plasmid or other non-pseudomolecule parts of assembly
        // println "${allowedTargetIdPattern} not matching $TNAME"
      } else if(markerMode && query_identity < 1) { //Not a 100% match, so for markers we check if no MM in last 3 bases - if notMarkerMode the required tag may not be present
          if(STRAND == '-' && (csTag =~ /^=[acgtnACGTN]{3,}/).count == 0){
            // println "minus 3'\t"+line
            kosher = false;
          } else if(STRAND == '+' && (csTag.substring(csTag.lastIndexOf('=')+1).matches('^[acgtnACGTN]{3,}\$')) == false) {
            // println "plus 3'\t"+line
            kosher = false;
          }
      }

      if(kosher) {
        def key = TNAME.replaceFirst("^(C|c)(H|h)(R|r)[_]?","")
        if(!scope.containsKey(key)) {
          scope << [(key) : []]
        }
        def feature = [
          "name" : QNAME,
          "value" : [ tstart+1, tend ],
          "evidence": [
            "identity" : query_identity,
            "coverage" : query_coverage,
            "strand": STRAND,
          ]
        ]
        if(csTag != null) {
          feature.evidence <<  ["cs" : csTag]
        }
        scope[key] << feature
      }
    }
}

//GROUP TOGETHER FEATURES FROM/IN SAME BLOCK
annotation.blocks = []
scope.each { k, features ->
  current = [ "scope": k, "featureType": "linear", "features": []]
  current.features = features
  annotation.blocks << current
}
//RECORD NUM FEATURS PER BLOCK
def total = 0;
counts.withWriterAppend{ wr ->
  annotation.blocks.each {
    wr.println annotation.name+"\t"+it.scope+"\t"+it.features.size()
    total += it.features.size()
  }
}
if(total == 0) {
  System.err.println('Zero sequences placed, terminating')
  System.exit(3)
}

if(output.endsWith('.gz')) {
  boolean append = false
  int BUFFER_SIZE = 8192
  Writer writer = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(out, append)), "UTF-8"), BUFFER_SIZE);
  try {
    writer.write(prettyPrint(toJson(annotation)))
    } catch (FileNotFoundException ex) {
      ex.printStackTrace();
    } catch (InterruptedException ex) {
      ex.printStackTrace();
    } catch (IOException ex) {
      ex.printStackTrace();
    } finally {
      try {
          writer.close();
      } catch (IOException ex) {
        ex.printStackTrace();
      }
    }
} else {
  out.text = prettyPrint(toJson(annotation))
}

if(outputTsv != null) {
  def outTsv = new File(outputTsv)
  outTsv.withWriterAppend{ wr ->
    wr.print "Name\tBlock\tstart\tend\tidentity\tcoverage\tstrand"
    wr.println csTag == null ? "" : "\tcs"
    // wr.println "Name\tBlock\tstart\tend\tidentity\tcoverage\tstrand\tcs"
    annotation.blocks.each { block ->
      block.features.each { feature ->
        cs = (feature.evidence.containsKey("cs") ? "\t"+feature.evidence.cs : "")
        wr.println feature.name+"\t"+block.scope+"\t"+feature.value[0]+"\t"+feature.value[1]+"\t"+feature.evidence.identity+"\t"+feature.evidence.coverage+"\t"+feature.evidence.strand+cs
      }
    }
  }
}

if(outputGff != null) {
  def outGff = new File(outputGff)
  outGff.withWriterAppend{ wr ->
    annotation.blocks.each { block ->
      block.features.each { feature ->
        gffEntry = [
          block.scope, basename, seqType, feature.value[0], feature.value[1],  ".", feature.evidence.strand, ".",
          "ID=${feature.name};Identity=${feature.evidence.identity};Coverage=${feature.evidence.coverage};"
        ].join('\t')
        wr.println gffEntry
      }
    }
  }
}
