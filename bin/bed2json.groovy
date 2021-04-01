#!/usr/bin/env groovy

import static groovy.json.JsonOutput.*
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream


@Grab('info.picocli:picocli:4.0.0-alpha-3') //command line interface
@Command(header = [
       //Font Name: Calvin S
       $/@|bold,blue  ╔╗ ╔═╗╔╦╗  ┌┬┐┌─┐  ╔═╗┬─┐┌─┐┌┬┐┌─┐┌─┐┬   |@/$,
       $/@|bold,blue  ╠╩╗║╣  ║║   │ │ │  ╠═╝├┬┘├┤  │ ┌─┘├┤ │   |@/$,
       $/@|bold,blue  ╚═╝╚═╝═╩╝   ┴ └─┘  ╩  ┴└─└─┘ ┴ └─┘└─┘┴─┘ |@/$
       ],
       description = "Convert BED to Pretzel JSON",
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
@Option(names = ["--make-private"], description = ["Make output dataset private (is public by default)"])
@Field private boolean makePrivate = false

@Option(names = ["--parent"], description = ["Parent dataset name"])
@Field private String parent

@Option(names = ["--name"], description = ["Basename used for dataset name"])
@Field private String name

// @Option(names = ["--namespace"], description = ["String used for dataset namespace"])
// @Field private String namespace

@Option(names = ["--short-name"], description = ["Short display name"])
@Field private String shortName

@Option(names = ["-b", "--bed"], description = ["Input BED file name"])
@Field private String bed = '/dev/stdin'

@Option(names = ["-a", "--allowed-target-id-pattern"], description = ["Provide target identifier patter if other than common chromosome naming"])
@Field private String allowedTargetIdPattern

@Option(names = ["-O", "--output"], description = ["JSON output file name"])
@Field private String output = '/dev/stdout'

@Option(names= ["-h", "--help"], usageHelp=true, description="Show this help message and exit.")
@Field private boolean helpRequested



//final int PAD4BASES = 9
//final int PAD4CHROMOSOMES = 5
final int BUFFER_SIZE = 8192
final String NEWLINE = System.lineSeparator();
final String SEP = '\t';

File bedFile = new File(bed)


def dataset = [:]
dataset.'public' = !makePrivate
dataset.name = name
dataset.namespace = "${parent}:${name}"
dataset.parent = parent

dataset.meta = [:]
// dataset.meta << ["source" : ""]
// dataset.meta << ["citation" : ""]
dataset.meta.shortName = shortName
dataset.blocks = []
//Process alignments
try {
  TreeMap scope = [:] //keep keys sorted as the corresponding blocks get displayed in order in pretzel

  def bedContentBuffer = new BufferedReader(new InputStreamReader(new FileInputStream(bedFile), "UTF-8"), BUFFER_SIZE);
  while ((line = bedContentBuffer.readLine()) != null && !line.isEmpty() && !line.matches('^@')) {
    def (BLOCK, START, END, ID, NUM, STRAND) = line.tokenize(',|\t')
    def position = STRAND == "+" ? [ START.toInteger() + 1, END.toInteger() ] : [ END.toInteger(), START.toInteger() + 1 ]
    key = BLOCK.replaceFirst("^(C|c)(H|h)(R|r)[_]?","")
      //
      //Skip non-chromosome blocks
      if(ID != "gap") {
        if(key.toLowerCase() =~ /^(chr|[0-9]|x|y|i|v)/ || (BLOCK =~ allowedTargetIdPattern) )  {
          if(!scope.containsKey(key)) {
            scope << [(key) : []]
          }
          scope[key] << ["name" : ID, "value" : position ] //, "strand": STRAND]
        }
      }
  }
    //GROUP TOGETHER FEATURES FROM/IN SAME BLOCK
    scope.each { k, features ->
      current = [ "scope": k, "featureType": "linear", "features": []]
      features.each { feature ->
        current.features << feature
      }
      dataset.blocks << current
    }
} catch (FileNotFoundException ex) {
  ex.printStackTrace();
} catch (InterruptedException ex) {
  ex.printStackTrace();
} catch (IOException ex) {
  ex.printStackTrace();
}

boolean append = false
if(output.endsWith('.gz')) {
  Writer  writer = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(output, append)), "UTF-8"), BUFFER_SIZE);
  try {
    writer.write(prettyPrint(toJson(dataset)))
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
    new File(output).withWriter { out ->
      out.println(prettyPrint(toJson(dataset)))
    }
}


