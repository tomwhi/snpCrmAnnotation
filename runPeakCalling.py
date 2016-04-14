#!/usr/bin/env python
# AUTHOR: Tom Whitington
# CREATE DATE: 8th July 2013

import cPickle, os, shutil, commands, pdb, sys
from optparse import OptionParser

import genomics.seqExpts.chipSeq as chipSeq
import genomics.formats.wigIO as wigIO
import genomics.seqExpts.core

def main():
    # Parse the command-line arguments...
    description = """usage: %prog [options] <chipAlignments> <controlAlignments>\n
Runs macs on the specified chip and control alignment (bam) files. Produces
a peak calls (bed) file, peak visualisation files, and an xml summary
file.

Inputs:
- ChIP-seq reads bam file.
- IGG control reads bam file. (optional - HACK)

Outputs:
- Writes all files to the specified output directory; "./peakCalling" by default
- Peaks bed file (a processed version of the macs peak files)
- Visualisation profiles:
-- Just for the peak regions and only for the ChIP sample.
-- Extended profile (wib, tdf)
-- Strand-specific profile (wib, tdf)
- XML summary of analysis (default stdout):
-- Metadata and statistics to report:
--- date
--- input files
--- output files
--- input paramters
--- number of peaks called
- All MACS output files
"""

    parser = OptionParser(usage = description)
    parser.add_option("--oc", dest = "outDir", default = "peakCalling",
                      help = "name of directory for output files will " + \
                          "not replace existing directory. [default: %default]")
    parser.add_option("--o", dest = "outDirNoReplace", default = None,
                      help = "name of directory for output files will " + \
                          "replace existing directory. [default: %default]")
    parser.add_option("--studyName", dest = "studyName", default = "UnknownStudy",
                      help = "Name of the study in which the experiment " + \
                          "took place. [default: %default]")
    parser.add_option("--chipChromatinName", dest = "chipChromatinName", default = None,
                      help = "Name of the chromatin tissue for the chip " + \
                          "sample. [default: %default]")
    parser.add_option("--controlChromatinName", dest = "controlChromatinName", default = None,
                      help = "Name of the chromatin tissue for the control " + \
                          "sample. [default: %default]")
    parser.add_option("--chipTargetName", dest = "chipTargetName", default = None,
                      help = "Name of the chromatin tissue for the chip " + \
                          "sample. [default: %default]")
    parser.add_option("--controlTargetName", dest = "controlTargetName", default = None,
                      help = "Name of the chromatin tissue for the control " + \
                          "sample. [default: %default]")
    parser.add_option("--chipCatalogNum", dest = "chipCatalogNum", default = None,
                      help = "Name of the chromatin tissue for the chip " + \
                          "sample. [default: %default]")
    parser.add_option("--controlCatalogNum", dest = "controlCatalogNum", default = None,
                      help = "Name of the chromatin tissue for the control " + \
                          "sample. [default: %default]")
    parser.add_option("--macsOptions", dest = "macsOptions",
                      default = "",
                      help = "Additional options for MACS. [default: %default]")
    parser.add_option("--genome", dest = "genome",
                      default = "hg19",
                      help = "The genome of the analysis. [default: %default]")
    parser.add_option("--singleProfExtn", dest = "singleProfExtn",
                      default = 1,
                      help = "Extension in bp for single-strand profiles. " + \
                          "[default: %default]")
    parser.add_option("--combinedProfExtn", dest = "combinedProfExtn",
                      default = 100,
                      help = "Extension in bp for combined profiles. " + \
                          "[default: %default]")
    parser.add_option("--baseDir", dest = "baseDir", default = None,
                      help = "The base directory of the chip-seq experiment." + \
                          " Setting this means that the chip experiment must" + \
                          " be located in this directory. It also means that" + \
                          " the analyses will be continued using the existing" + \
                          " analyser if it exists, and the analyser will be " + \
                          " \"pickled\" to the base directory afterwards.")
    parser.add_option("--debug", action="store_true", dest="debug",
                      help = "Debug the program using pdb.")
    (options, args) = parser.parse_args()

    # Parse the input parameters...

    if (options.debug):
        pdb.set_trace()

    if (len(args) < 1):
        print >> sys.stderr, "WRONG # ARGS: ", len(args)
        parser.print_help()
        sys.exit(1)

    chipAlnFilename = args[0]
    if len(args) == 2:
        controlAlnFilename = args[1]
    else:
        controlAlnFilename = None

    # Determine the output directory; create/over-write if needed...
    outDirName = None
    if options.outDirNoReplace != None:
        # Told not to replace the named output folder if it exists.

        outDirName = options.outDirNoReplace
        if os.direxists(os.path.join(os.getcwd()), outDirName):
            # Quit, since the folder already exists:
            print >> sys.stderr, "Output folder already exists:", outDirName
            sys.exit(1)
    else:
        outDirName = options.outDir

    if options.baseDir != None:
        # Try to determine the chip target and chromatin values from the
        # base directory absolute location:
        # FIXME: This will not work if the specified base directory is
        # not part of my planned ChIP-seq directory heirarchy.
        chipBamAbsLoc = os.path.abspath(chipAlnFilename)
        chipBamToks = chipBamAbsLoc.split("/")
        chipChromatinName = chipBamToks[-4]
        chipTargetName = chipBamToks[-3]
        if controlAlnFilename != None:
            controlBamAbsLoc = os.path.abspath(controlAlnFilename)
            controlBamToks = controlBamAbsLoc.split("/")
            controlChromatinName = controlBamToks[-4]
            controlTargetName = controlBamToks[-3]
        else:
            controlChromatinName = None
            controlTargetName = None
    else:
        # Use the specified values:
        chipChromatinName = options.chipChromatinName
        controlChromatinName = options.controlChromatinName
        chipTargetName = options.chipTargetName
        controlTargetName = options.controlTargetName

    # Objects representing the chip and control protocols that were carried out:
    chipChromatin = chipSeq.Chromatin(chipChromatinName)
    chipAntibody = chipSeq.Antibody(chipTargetName, options.chipCatalogNum)
    chipProtocol = chipSeq.Chip(chipChromatin, chipAntibody)
    chipLibrary = genomics.seqExpts.core.Library(chipProtocol)

    if controlAlnFilename != None:
        controlChromatin = chipSeq.Chromatin(controlChromatinName)
        controlAntibody = chipSeq.Antibody(controlTargetName, options.controlCatalogNum)
        controlProtocol = chipSeq.Chip(controlChromatin, controlAntibody)
        controlLibrary = genomics.seqExpts.core.Library(controlProtocol)
    else:
        controlLibrary = None

    # The experiment:
    expt = chipSeq.ChipseqExpt(chipLibrary, controlLibrary,
                               studyName=options.studyName)

    # Set the two input bam files for the experiment:
    expt.setChipAln(chipAlnFilename)
    expt.setControlAln(controlAlnFilename)

    # Set up the combinations of parameters for generating profiles:
    paramsForProfs = []
    singleProfExtn = int(options.singleProfExtn)
    combinedProfExtn = int(options.combinedProfExtn)
    if singleProfExtn > 0:
        plusStrandProf = wigIO.ProfileParams(singleProfExtn,
                                             wigIO.ProfileParams.plus)
        minusStrandProf = wigIO.ProfileParams(singleProfExtn,
                                              wigIO.ProfileParams.minus)
        paramsForProfs.append(plusStrandProf)
        paramsForProfs.append(minusStrandProf)
    if combinedProfExtn > 0:
        combinedStrandProf = wigIO.ProfileParams(combinedProfExtn,
                                                 wigIO.ProfileParams.combined)
        subtractedStrandProf = wigIO.ProfileParams(combinedProfExtn,
                                                   wigIO.ProfileParams.subtracted)
        paramsForProfs.append(combinedStrandProf)
        paramsForProfs.append(subtractedStrandProf)

    # The peak calling analysis:
    analysis = chipSeq.PeakCallingAna(expt, "PeakCalling", outDirName,
                                      paramsForProfs, genome=options.genome,
                                      macsOptionsStr=options.macsOptions)

    # Create the analyser, for running the analysis. Consider
    # the base directory if told to do so:
    analyser = None

    pickledAnalyserFilename = None
    if options.baseDir != None:
        # Continue from existing analyser if it's there:
        # FIXME: Hard-coded filename:
        pickledAnalyserFilename = options.baseDir + "/Analyser.pk"
        pickledAnalyserFile = None
        try:
            pickledAnalyserFile = open(pickledAnalyserFilename)
        except IOError:
            pass

        # Check all the paths to make sure they're consistent with
        # the specified base directory...

        # FIXME: The following directory-checking rules are quite
        # in-depth and don't seem to belong here. Should they be inside
        # the anchoredSeqExtpAnalyser class? I don't think so, since it should
        # not have to know about all the various component analyses.

        # FIXME: This code currently only works if the input files, output
        # directory and base directory are all specified as absolute paths. This
        # seems restrictive and buggy.

        # FIXME: Consider relative paths in all the following
        # file name comparisons...
        if pickledAnalyserFile != None:
            analyser = cPickle.load(pickledAnalyserFile)
            if (analyser.baseDir != options.baseDir):
                print >> sys.stderr, "Previous analyser's base directory " + \
                    "\"" + analyser.getBaseDir() + "\" is not the same as " + \
                    "the specified base directory \"" + options.baseDir + "\"."
                sys.exit(1)
        else:
            analyser = genomics.seqExpts.core.AnchoredSeqExptAnalyser(options.baseDir)

        # Check that the chip bam file is located inside the base dir:
        lastSlashIdx = chipAlnFilename.rfind("/")
        chipAlnDir = chipAlnFilename[:lastSlashIdx]

        if (chipAlnDir != options.baseDir):
            print >> sys.stderr, "Directory of ChIP bam file " + \
                "\"" + chipAlnDir + "\" is not the same as " + \
                "the specified base directory \"" + options.baseDir + "\"."
            sys.exit(1)

        # Check that the output directory is inside the base dir:
        lastSlashIdx = outDirName.rfind("/")
        outDirBase = outDirName[:lastSlashIdx]
        if (outDirBase != options.baseDir):
            print >> sys.stderr, "Directory of output dir " + \
                "\"" + outDirName + "\" is not the same as " + \
                "the specified base directory \"" + options.baseDir + "\"."
            sys.exit(1)
    else:
        # The analyser does not have a base directory for the ChIP-seq
        # experiment:
        analyser = genomics.seqExpts.core.SeqExptAnalyser()

    # Delete the output directory if it exists, and re-make it:
    shutil.rmtree(outDirName, ignore_errors=True)
    os.makedirs(outDirName)

    analyser.addAnalysis(analysis)

    # Run the analysis:
    analyser.runAnalyses()

    # Generate XML:
    analyser.writeXML(sys.stdout)

    if options.baseDir:
        # Pickle the analyser to carry out future analyses. This
        # will over-write the previous pickling of the analyser, if
        # it existed:
        pickledAnalyserFile = open(pickledAnalyserFilename, 'w')
        cPickle.dump(analyser, pickledAnalyserFile)
        pickledAnalyserFile.close()


if __name__ == "__main__":
    sys.exit(main())
