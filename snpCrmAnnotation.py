#!/usr/bin/env python

import commands, MySQLdb, pdb, sys, cPickle
import getSeqFa, snvAnnotation
from optparse import OptionParser


def motEffectCmp(effect1, effect2):
    if effect1.absDiff < effect2.absDiff:
        return -1
    elif effect1.absDiff == effect2.absDiff:
        return 0
    else:
        return 1


def main():
    # Parse the command-line arguments...
    description = """usage: %prog [options] <snvsVcfFile> <refGenomeFasta> <outFile>\n
Inputs:
- VCF file of SNVs to annotate (- specifies stdin)
- Access to database containing LD information and all feature annotations.
- Parameters specifying which features to consider for overlap.

Outputs:
- CHANGED: Now just outputs a pickled tuple containing the snv annotation info.
- SNVs that had one or more SNVs in LD overlapped by the specified combination of features.
- Definition of "overlap":
-- Need to allow any combination of these factors:
--- Broad region specified by tissue name
--- Footprint specified by tissue name
--- Affected motif
---- Optionally reporting phylogenetic conservation of the best motif hit
---- - AND/OR -
---- Optionally reporting corresponding ChIP-seq peak overlap, for those motifs that correspond to a TF analysed with ChIP-seq
--- ChIP-seq peak, with flexible summit distance cutoff
--- Phylogenetic conservation of the SNV itself.
"""

    lowCons = -10 # Lower than lowest possible conservation score; marker.

    parser = OptionParser(usage = description)
    parser.add_option("--mysqlUname", dest = "uname", default = "tom",
                      help = "Mysql username for logging into the snvEffect" + 
                      "database. [default: %default]")
    parser.add_option("--chipseqBaseDir", dest = "chipseqBaseDir",
                      default = "/mnt/skjoldNoBackup/pca_brca/",
                      help = "Base directory for ChIP-seq. [default: %default]")
    parser.add_option("--motifBaseDir", dest = "motifBaseDir",
                      default = "/home/tom/Work/data/PWMs/db/",
                      help = "Base directory for PWMs. [default: %default]")
    parser.add_option("--mysqlHost", dest = "host", default = "localhost",
                      help = "Mysql host server to use. [default: %default]")
    parser.add_option("--r2thresh", dest = "r2thresh", default = "0.3",
                      help = "R-squared threshold for determining" + 
                      " LD. [default: %default]")
    parser.add_option("--distThresh", dest = "distThresh", default = "-1",
                      help = "Distance threshold for determining" + 
                      " associated SNVs. Use a positive value to " + 
                      "override r2 threshold. [default: %default]")
    parser.add_option("--mysqlSelDb", dest = "snvEffectDb", default = "snvEffect",
                      help = "Mysql database name for logging into the snvEffect" + 
                      "database.")
    parser.add_option("--trace", action = "store_true", dest = "trace",
                      help = "Print out detailed info on overlaps.")
    parser.add_option("--regRegTissue", dest = "regRegTissue", default = None,
                      help = "Tissue for regulatory region overlap." +
                      " [default: %default]")
    parser.add_option("--footprintTissue", dest = "footprintTissue",
                      default = None,
                      help = "Tissue for DNase HS footprint overlap." +
                      " [default: %default]")
    parser.add_option("--checkMotif", dest="checkMotif", default = -1,
                      help = "Check for overlap of affected motifs at " +
                      "the SNV. Use the specified threshold to determine if " +
                      "a given motif hit is affected." +
                      " [default: %default]")
    parser.add_option("--outputPvals", action="store_true", dest="outputPvals",
                      help = "Output motif scores as p-values and pval ratios.")
    parser.add_option("--showBestHitSeq", action = "store_true",
                      dest = "showBestHitSeq",
                      help = "Show sequence of best allele in motif" + 
                      " hits. NOTE: Not really needed by program, but needed " +
                      "by scoreSNPs module interface. Should fix this.")
    parser.add_option("--motifMinScoreDiff", dest="scoreDiffThresh", default = 0.5,
                      help = "Minimum difference in log10 best p-values for " +
                      "best- vs worst-scoring alleles. [default: %default]")
    parser.add_option("--motifMinBestScore", dest="bestScoreThresh", default = 3,
                      help = "Minimum score of best allele for snv motif effects." +
                      " [default: %default]")
    parser.add_option("--maxMotifOverlaps", dest="maxMotifOverlaps",
                      default = -1,
                      help = "Maximum number of motif overlaps to report " +
                      "for each SNV. [default: %default]")
    parser.add_option("--motifCons", dest="motifCons", default = lowCons,
                      help = "Check for conservation of affected motifs." +
                      " Use the specified threshold to determine conservation."+
                      " [default: %default]")
    parser.add_option("--matchPeak", dest="matchPeak", default = -1,
                      help = "Check if affected motif hits are occupied." + 
                      " Use the specified peak distance threshold." +
                      " [default: %default]")
    parser.add_option("--onlyOccupied", action = "store_true", dest = "onlyOccupied",
                      help = "Only report occupied affected motif hits.")
    parser.add_option("--snvCons", dest="snvCons", default = lowCons,
                      help = "Check if the SNV itself is conserved." +
                      " Use the specified threshold to determine conservation."+
                      " [default: %default]")
    parser.add_option("--motifPrefixes", dest="motifPrefixes",
                      default = "/home/tom/Work/data/PWMs/allMotifPrefixes.txt",
                      help = "List of prefixes for the motifs to be analysed." +
                      " [default: %default]")
    parser.add_option("--motifPrefixes_DNase", dest="motifPrefixes_DNase",
                      default = None,
                      help = "Prefixes for the motifs to be analysed at DNase-." +
                      "overlapped SNPs. [default: %default]")
    parser.add_option("--dnaseFocus", dest="dnaseFocus",
                      default = None,
                      help = "Name of DNase experiment that is particularly." +
                      "relevant to this annotation. Used in PWM scanning. [default: %default]")
    parser.add_option("--bigWigs", dest="bigWigs",
                      default = None,
                      help = "BigWig files for SNP annotation, comma-separated list of." +
                      "<variableName>=<filename> pairs. [default: %default]")
    parser.add_option("--peaksTabixFile", dest="peaksTabixFile",
                      default = None,
                      help = "Tabix file containing chip-seq peaks locations." +
                      " [default: %default]")
    parser.add_option("--bedFeatures", dest="bedFeatures",
                      default = "/home/tom/Work/data/hg19/ucsc/knownGene_ExonLocations_hg19.bed",
                      help = "A comma-separated list of bed files designating additional " +
                      "features to overlap the input SNPs." +
                      " [default: %default]")
    parser.add_option("--chipseqTabixFile", dest="chipseqTabixFile",
                      default = "/home/tom/Work/data/hg19/InternalChipseq/chipseqBindPosTable_hg19.bed.gz",
                      help = "Tabix file containing estimated chip-seq bind pos locations." +
                      " [default: %default]")
    parser.add_option("--consDir", dest="consDir", default = None,
                      help = "Wib directory containing conservation data.")
    parser.add_option("--peak", dest="peak", default = -1,
                      help = "Look in the specified range for ChIP-seq peaks." +
                      " Use the specified peak distance threshold." +
                      " [default: %default]")
    parser.add_option("--tmpDir", dest="tmpDir", default = "/tmp",
                      help = "Location of the temporary directory, where tmp files will be stored." +
                      " [default: %default]")
    parser.add_option("--chipseqDist", dest="chipseqDist", default = -1,
                      help = "Look in the specified range for subtracted profile." +
                      " Binding positions. Use the specified peak distance threshold." +
                      " -1 => Don't use. Only used if peaks are too. [default: %default]")
    parser.add_option("--debug", action="store_true", dest="debug",
                      help = "Debug the program using pdb.")
    (options, args) = parser.parse_args()

    # Parse the input parameters...

    if (options.debug):
        pdb.set_trace()

    # Make sure the 1 required input file exists:
    if (len(args) != 3):
        print >> sys.stderr, "WRONG # ARGS: ", len(args)
        parser.print_help()
        sys.exit(1)

    # Get all SNVs to annotate. Will process all those specified on stdin:
    snvsFilename = args[0]

    # Get the sequence retriever:
    seqGetter = getSeqFa.seqRetriever(args[1])

    # Annotate the SNVs:
    snvAnnotator = snvAnnotation.SNV_Annotator(seqGetter, options)
    annotGroup = snvAnnotator.annotateSNVs(snvsFilename)

    # Output annotation information, in a parseable form that will
    # make it easy to generate visualisations...
    # Currently just pickling the annotation object:
    outfile = open(args[2], 'w')
    annotGroup.pickle(outfile)
    outfile.close()


if __name__ == "__main__":
    sys.exit(main())
