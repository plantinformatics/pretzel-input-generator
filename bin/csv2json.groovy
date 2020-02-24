#!/usr/bin/env groovy

import static groovy.json.JsonOutput.*
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream


@Grab('info.picocli:picocli:4.0.0-alpha-3') //command line interface
@Command(header = [
       //Font Name: Calvin S
       $/@|bold,blue  ╔═╗╔═╗╦  ╦  ┌┬┐┌─┐  ╔═╗┬─┐┌─┐┌┬┐┌─┐┌─┐┬   |@/$,
       $/@|bold,blue  ║  ╚═╗╚╗╔╝   │ │ │  ╠═╝├┬┘├┤  │ ┌─┘├┤ │   |@/$,
       $/@|bold,blue  ╚═╝╚═╝ ╚╝    ┴ └─┘  ╩  ┴└─└─┘ ┴ └─┘└─┘┴─┘ |@/$
       ],
       description = "Convert CSV genetic map to Pretzel dataset JSON",
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
@Option(names = ["-p", "--private"], description = ["Make output dataset private (is public by default)"])
@Field private boolean makePrivate = false

@Option(names = ["-c", "--csv"], description = ["Input CSV file name"])
@Field private String csv = '/dev/stdin'

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

File csvFile = new File(csv)

//Process alignments
try {
  def dataset = [:]
  dataset.public = !makePrivate
  dataset.meta = [:]
  dataset.meta << ["shortName" : ""]
  dataset.meta << ["source" : ""]
  dataset.meta << ["citation" : ""]

  dataset.name = ""
  dataset.namespace = ""
  dataset.parent = ""
  dataset.blocks = []
  TreeMap scope = [:] //keep keys sorted as the corresponding blocks get displayed in order in pretzel

  def csvContentBuffer = new BufferedReader(new InputStreamReader(new FileInputStream(csvFile), "UTF-8"), BUFFER_SIZE);
  while ((line = csvContentBuffer.readLine()) != null && !line.isEmpty() && !line.matches('^@')) {
    def (MARKER, BLOCK, POS1, POS2) = line.tokenize(',|\t')
    def POS1NUM = POS1.isInteger() ? POS1.toInteger() : POS1.toDouble()
    def POS2NUM = POS2.isInteger() ? POS2.toInteger() : POS2.toDouble()
    def position = POS2 == null ? [POS1NUM] : [POS1NUM, POS2NUM]


    key = BLOCK.replaceFirst("^(C|c)(H|h)(R|r)[_]?","")
      //Skip non-chromosome blocks
      if(key.toLowerCase() =~ /^(chr|[0-9]|x|y|i|v)/ || (BLOCK =~ allowedTargetIdPattern) )  {
        if(!scope.containsKey(key)) {
          scope << [(key) : []]
        }
        scope[key] << ["name" : MARKER, "value" : position]
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
  new File(output).withWriter { out ->
    out.println(prettyPrint(toJson(dataset)))
  }
} catch (FileNotFoundException ex) {
  ex.printStackTrace();
} catch (InterruptedException ex) {
  ex.printStackTrace();
} catch (IOException ex) {
  ex.printStackTrace();
}



