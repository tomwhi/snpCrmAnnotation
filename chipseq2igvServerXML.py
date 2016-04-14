#!/usr/bin/env python

import cPickle, datetime, os, pdb, re, sys
from optparse import OptionParser
import xml.etree.ElementTree as ET

import genomics.seqExpts.chipSeq as chipSeq
import genomics.formats.wigIO as wigIO
import genomics.seqExpts.core


def generateAnalyserXML(parentNode, analysers, baseDataDir2urlBase):
    """Generates XML content for the specified IGV data registry xml document,
    for the analysers (i.e. not for specific actuall analyses)."""

    # Make a node for the bam files:
    bamFilesNode = ET.SubElement(parentNode, "Category")
    bamFilesNode.set("name", "BamFiles")

    # Generate XML content for the individual analysers...
    for analyser in analysers:
        # FIXME:
        # Implementing here rather than in analyser class, which is a bit
        # naughty. Might need to re-factor to fix...

        # Get the chip bam file location:
        analysis = analyser.getAnalyses().values()[0]
        absChipBamLoc = os.path.abspath(analysis.getExpt().getChipAln())
        urlLoc = absChipBamLoc.replace(baseDataDir2urlBase[0],
                                       baseDataDir2urlBase[1])

        currAnalysisNode = ET.SubElement(bamFilesNode, "Resource")

        currAnalysisNode.set("name", analysis.getExpt().getName())
        currAnalysisNode.set("path", urlLoc)


def generatePeakCallingXML(parentNode, analysers, baseDataDir2urlBase):
    """Generates XML content for the specified IGV data registry xml document,
    for the peak-calling analyses."""

    # Get the peak-calling analyses:
    peakCallingAnalyses = \
        map(lambda analyser: analyser.getAnalysis("PeakCalling"), analysers)

    # Generate the base-level nodes of the XML document for the
    # peak calling...
    allPeakCallingsNode = ET.SubElement(parentNode, "Category")
    allPeakCallingsNode.set("name", "PeakCallingAnalyses")

    allPeaksNode = ET.SubElement(allPeakCallingsNode, "Category")
    allPeaksNode.set("name", "Peaks")

    allProfilesNode = ET.SubElement(allPeakCallingsNode, "Category")
    allProfilesNode.set("name", "Profiles")

    # Currently, assume all specified chip-seq peak-calling analyses have the
    # same profile parameters. Hence, retrieve the profile parameters from
    # (arbitrarily) the first chip-seq peak-calling analysis:
    profParams = peakCallingAnalyses[0].getProfileParams()

    profParams2paramsXmlNodes = {}

    for params in profParams:
        paramsStr = params.toString()
        currParamsNode = ET.SubElement(allProfilesNode, "Category")
        currParamsNode.set("name", paramsStr)
        profParams2paramsXmlNodes[params] = currParamsNode

    # Sort the analyses so that the tracks will appear in a sensible default
    # order when loaded...
    peakCallingAnalyses.sort(chipSeq.cmpAnalyses)

    # Generate XML content for the individual peak calling analyses, editing
    # the XML skeleton established above:
    for analysis in peakCallingAnalyses:
        exptName=analysis.getExpt().getName()

        # Get the peak-calls base XML node and add XML to it for this analysis:
        analysis.writePeakCallsIGV_XML(allPeaksNode,
                                       baseDataDir2urlBase, exptName)

        # Generate XML for all the profiles for this analysis:
        for params in profParams:
            # NOTE: A precondition to this function is that all analysers'
            # chipseqs have the same profile parameters. So, I can assume
            # that I'll find a profile given the current params, in the current
            # analysis:
            currProfile = analysis.getProfile(params)

            # Get the base XML profile for these parameters, and add XML to it
            # by calling the Profile's method. Tell it a tup showing how to
            # convert from base dir to url:
            currParamsXMLnode = profParams2paramsXmlNodes[params]

            # Set the starting track attributes (colour, default data limits
            # etc.) based on the profile parameters and the type of ChIP-seq
            # experiment:
            
            # FIXME: This is another part that should be re-factored when I'm
            # figuring out how to generate XML properly. For the time being,
            # I will call a function here, to construct a simple
            # attribute->value dictionary recording the various attributes for
            # the given profile:
            att2val = getAtt2Val(analysis, params)

            currProfile.writeIGV_XML(currParamsXMLnode, baseDataDir2urlBase,
                                     exptName, att2val)


def getAtt2Val(chipseqAnalysis, params):
    """Temporary/yucky function (but then again so is this whole program).
    Returns an attribute to string value mapping for a given wiggle profile
    track, based on the analysis features (such as chromatin and name) and
    parameters.

    Rules determining the default colour for a given ChIP-seq track:
    Histone modification (H[1234].*) -> dark purple
    Polymerase or mediator (POLR2A or MED[0-9]+) -> bright green
    CTCF/RAD21/SMC[0-9]+.+ -> orange
    Otherwise, assume it's a transcription factor -> pink
    """

    att2val = {}

    currChip = chipseqAnalysis.getExpt().getLibSeq(chipSeq.ChipseqExpt.chipName).getLibrary().getProtocol()
    targetStr = currChip.antibody.getTargetName()

    # Implementing the colour rules described above
    colour = "255,115,255"
    if re.match("H[1234].*", targetStr) != None:
        colour = "133,7,128"
    elif re.match("MED[0-9]+", targetStr) != None or re.match("POLR2A", targetStr) != None:
        colour = "74,255,77"
    elif re.match("CTCF", targetStr) != None or re.match("RAD21", targetStr) != None or \
            re.match("SMC[0-9]+.*", targetStr) != None:
        colour = "255,232,79"

    # Always set profile tracks to use max as the default windowing function:
    att2val["windowFunction"] = "max"

    att2val["color"] = colour

    # Always use default height of 20:
    att2val["height"] = "20"

    att2val["trackLine"] = "viewLimits=0:50"
    return att2val


def main():
    # Parse the command-line arguments...
    description = """usage: %prog [options] <chipseqAnalysers>\n
Inputs:
A list of ChIP-seq analysers, comma-separated. Each file location is a path to
a pickled analyser object.

Outputs:
XML to standard output, which can then be loaded by an IGV server as a data
registry file.
"""

    parser = OptionParser(usage = description)
    parser.add_option("--debug", action="store_true", dest="debug",
                      help = "Debug the program using pdb.")
    (options, args) = parser.parse_args()

    # Parse the input parameters...

    if (options.debug):
        pdb.set_trace()

    if (len(args) != 1):
        print >> sys.stderr, "WRONG # ARGS: ", len(args)
        parser.print_help()
        sys.exit(1)

    chipseqAnalyserLocs = args[0].strip(",").split(",")

    # Parse the chip-seq experiment analysis pk files into analyser objects:
    analysers = map(lambda analyserLoc: cPickle.load(open(analyserLoc)), chipseqAnalyserLocs)

    # Produce a root node for the IGV data registry XML document:
    xmlRegistryNode = ET.Element('Global')
    xmlRegistryNode.set("name", "ChIP-seq data - " + str(datetime.datetime.now()))
    xmlRegistryNode.set("infolink", "http://localhost")
    xmlTree = ET.ElementTree(xmlRegistryNode)

    # FIXME: Hacky and bug-prone. A tuple for converting base file location
    # into the URL base:
    baseDataDir2urlBase = ("/home/tom/Work/data/processed/hg19/chipseq",
                           "http://localhost/chipseqData")

    generateAnalyserXML(xmlRegistryNode, analysers, baseDataDir2urlBase)

    # Produce xml for the peak-calling analyses:
    generatePeakCallingXML(xmlRegistryNode, analysers, baseDataDir2urlBase)

    # Write the XML document to sys.stdout:
    xmlTree.write(sys.stdout)


if __name__ == "__main__":
    sys.exit(main())
