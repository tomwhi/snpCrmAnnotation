#!/usr/bin/env python

import os, math, numpy, pdb, re, sys
import genomics.formats.bed as bed
import xml.etree.ElementTree as ET
import genomics.seqExpts.core as core
import genomics.formats.conversion as conversion
import genomics.formats.wigIO as wigIO
import matplotlib.pyplot as plt
import cmdlineProgs, utility


# FIXME: Later, perhaps implement an option to parse this information from
# XML, for the Chip, Chromatin, and Antibody classes.

class Chip(core.Protocol):
    """A single chromatin immunoprecipitation."""

    def __init__(self, chromatin, antibody):
        core.Protocol.__init__(self)

        self.chromatin = chromatin
        self.antibody = antibody

    def getName(self):
        # Currently just concatenate chromatin tissue and antibody target:
        return self.chromatin.getTissueName() + "_" + \
            self.antibody.getTargetName()


class Chromatin:
    """A single chromatin preparation."""

    # FIXME: Later, perhaps implement an option to parse this information from
    # XML.

    def __init__(self, tissueName):
        self.tissueName = tissueName

    def getTissueName(self):
        if self.tissueName == None:
            return "UnknownTissue"
        else:
            return self.tissueName


class Antibody:
    """A single antibody."""

    def __init__(self, targetName, catalogNum):
        self.targetName = targetName
        self.catalogNum = catalogNum

    def getTargetName(self):
        if self.targetName == None:
            return "UnknownTarget"
        else:
            return self.targetName


class ChipseqExpt(core.SequencingExpt):
    """A single ChIP-seq experiment."""

    chipName = "chip"
    controlName = "control"
    
    def __init__(self, chipLibrary, controlLibrary, studyName="UnknownStudy"):
        core.SequencingExpt.__init__(self)

        # Create empty LibrarySequencing objects representing the chip and control
        # sequence sets:
        chip_LibSeq = core.LibrarySequencing(chipLibrary)
        control_LibSeq = core.LibrarySequencing(controlLibrary)
        self.setLibSeq(ChipseqExpt.chipName, chip_LibSeq)
        self.setLibSeq(ChipseqExpt.controlName, control_LibSeq)
        self.studyName = studyName

    def setChipAln(self, chipAlnFilename):
        self.getLibSeq(ChipseqExpt.chipName).setBam(chipAlnFilename)

    def setControlAln(self, controlAlnFilename):
        self.getLibSeq(ChipseqExpt.controlName).setBam(controlAlnFilename)

    def getChipAln(self):
        """Returns a string representing the location of the bam file of chip
        read alignments for this experiment."""
        return self.getLibSeq(ChipseqExpt.chipName).getBam()

    def getControlAln(self):
        """Returns a string representing the location of the bam file of chip
        read alignments for this experiment."""
        return self.getLibSeq(ChipseqExpt.controlName).getBam()

    def getName(self):
        """Returns a string for a name for this experiment."""

        # FIXME: Need to refactor when I've thought about how to do this
        # properly. At the moment, just generate something from
        # the target, chromatin, and study:

        # FIXME: I want to have the experiment accession shown in the experiment
        # name too. However, (I think) I have not stored it in the analysis
        # or experiment object. So, I will hack it out here, by assuming it's
        # in the path of the chip filename:
        chipBamLoc=self.getLibSeq(ChipseqExpt.chipName).getBam()
        exptAcc=chipBamLoc.split("/")[-2]

        return self.studyName + "_" + \
            self.getLibSeq(ChipseqExpt.chipName).getLibrary().getProtocol().getName() + \
            "_" + exptAcc


# FIXME: Move this inside class Antibody before re-running all peak-calling:
def getTargetInt(targetStr):
    # FIXME: Nasty and hard-coded determination of the "type" of an antibody
    # target.
    if re.match("H[1234].*", targetStr) != None:
        return 1
    elif re.match("MED[0-9]+", targetStr) != None or re.match("POLR2A", targetStr) != None:
        return 2
    elif re.match("CTCF", targetStr) != None or re.match("RAD21", targetStr) != None or \
            re.match("SMC[0-9]+.*", targetStr) != None:
        return 3
    else:
        # It must be a TF.
        return 0
    

def cmpAnalyses(ana1, ana2):
    """Compares two chip-seq analyses. Rules for comparison:

    Consider target in terms of it's class (histone, polII, TF, etc.),
    and then in terms of the lexicographical comparison of the actual
    target string.

    Then consider the following attributes lexicographically:
    - Tissue
    - Study
    - Experiment
"""
    ana1expt = ana1.getExpt()
    ana2expt = ana2.getExpt()

    ana1exptAcc = \
        ana1expt.getLibSeq(ChipseqExpt.chipName).getBam().split("/")[-2]
    ana1target = \
        ana1expt.getLibSeq(ChipseqExpt.chipName).getLibrary().getProtocol().antibody.getTargetName()
    ana1tissue = \
        ana1expt.getLibSeq(ChipseqExpt.chipName).getLibrary().getProtocol().chromatin.getTissueName()
    ana1study = ana1expt.studyName

    ana2exptAcc = \
        ana2expt.getLibSeq(ChipseqExpt.chipName).getBam().split("/")[-2]
    ana2target = \
        ana2expt.getLibSeq(ChipseqExpt.chipName).getLibrary().getProtocol().antibody.getTargetName()
    ana2tissue = \
        ana2expt.getLibSeq(ChipseqExpt.chipName).getLibrary().getProtocol().chromatin.getTissueName()
    ana2study = ana2expt.studyName
        
    if getTargetInt(ana1target) < getTargetInt(ana2target):
        return -1
    elif getTargetInt(ana1target) > getTargetInt(ana2target):
        return 1
    else:
        if ana1target < ana2target:
            return -1
        elif ana1target > ana2target:
            return 1
        else:
            if ana1tissue < ana2tissue:
                return -1
            elif ana1tissue > ana2tissue:
                return 1
            else:
                if ana1study < ana2study:
                    return -1
                elif ana1study > ana2study:
                    return 1
                else:
                    if ana1exptAcc < ana2exptAcc:
                        return -1
                    elif ana1exptAcc > ana2exptAcc:
                        return 1
                    else:
                        return 0


class ChipseqAnalysis(core.SeqAnalysis):
    """An (abstract) chip-seq analysis."""

    def __init__(self, expt, name, outDir):
        # Precondition: expt must be a ChipseqExpt instance

        core.SeqAnalysis.__init__(self, expt, name, outDir)

    def editIO_XML(self, xmlElem):
        """Generate XML content specific to a ChIP-seq analysis."""

        core.SeqAnalysis.editIO_XML(self, xmlElem)

        # Experiment? Bam file? Or what?
        # FIXME: Not sure what content (if anything) I should add here.

    def editParams_XML(self, xmlElem):
        """Generate XML parameters content specific to a ChIP-seq
        analysis."""

        core.SeqAnalysis.editParams_XML(self, xmlElem)


class PeakCallingAna(ChipseqAnalysis):
    """A ChIP-seq peak calling analysis."""

    def __init__(self, expt, name, outDir, paramsForProfs, genome="hg19",
                 macsOptionsStr = ""):
        ChipseqAnalysis.__init__(self, expt, name, outDir)

        self.genome = genome

        self.macsOptionsStr = macsOptionsStr

        # The peaks bed file, produced after running MACS and converting from
        # MACS' xls file:
        self.peaksBedFilename = None

        # paramsForProfs is an array of zero or more ProfileParams objects,
        # indicating what kind of profiles to generate.

        # This object will be used to generate all profiles for
        # peak visualisations:
        profPrefix = outDir + "/Profile"

        # Name the profile according to the Chip LibrarySequencing
        chipLibSeq = expt.getLibSeq(ChipseqExpt.chipName)
        profName = chipLibSeq.getName()
        self.profileGen = wigIO.ProfileGenerator(name, profPrefix,
                                                 paramsForProfs,
                                                 genome)

    def getProfileParams(self):
        return self.profileGen.getAllParams()

    def getProfile(self, params):
        return self.profileGen.getProf(params.getProfType(),
                                       params.getExtension())

    def run(self):
        """Assumes that the chip and control mapping attributes have been set
        in the analysed experiment. Runs the analysis."""

        ChipseqAnalysis.run(self)

        self.runMACS()
        self.generateProfiles()

    def runMACS(self):
        """Runs MACS with the specified parameters."""

        # Retrieve the chip and control mapping file locations:
        chipBam = self.expt.getChipAln()
        controlBam = self.expt.getControlAln()

        # Run MACS:
        chipLibSeq = self.expt.getLibSeq(ChipseqExpt.chipName)
        chipLibName = chipLibSeq.getName()
        cmdlineProgs.run_MACS(chipBam, controlBam, self.outDir,
                              name=chipLibName,
                              macsOptionsStr = self.macsOptionsStr,
                              verbose=True)

        # Check that the peaks.xls file exists and quit, reporting an error,
        # if it does not:
        macsPeaksXlsFilename = self.outDir + "/" + chipLibName + "_peaks.xls"
        try:
            with open(macsPeaksXlsFilename): pass
        except IOError:
            print 'Error: MACS failed to produce a peaks.xls file. Quitting.'
            sys.exit(1)

        # Convert the MACS output to bed format and save the name of the
        # resulting file:
        # FIXME: Presently hard-coded:
        # FIXME: Are these paths correct, given MACS' interface?:
        self.peaksBedFilename = self.outDir + "/peaks.bed"
        conversion.macs2bed(macsPeaksXlsFilename, self.peaksBedFilename)

    def generateProfiles(self):
        """Generates the profiles for visualising the peaks, and convert them
        to the desired binary formats."""

        # Generate all profiles for all regions:
        peaksBedTrack = bed.BED_Track(open(self.peaksBedFilename))
        self.profileGen.makeAllRegProfiles(self.getExpt().getChipAln(),
                                           peaksBedTrack)

        # The profile generator can be retrieved later in order to get the
        # locations of the resulting wig/wib/tdf files.

    def editIO_XML(self, xmlElem):
        """Generate XML content specific to a peak-calling analysis."""

        ChipseqAnalysis.editIO_XML(self, xmlElem)
        
        # FIXME: Write chip and control bam file nodes, peaks outfile node, and outfile nodes for each of the wib and tdf profile files produced. NOTE: This would be a good place to automatically and cleanly extract information from when attempting to produce an xml file for the IGV server.
        peakCallingNode = ET.SubElement(xmlElem, 'PeakCallingIO')

    def editParams_XML(self, xmlElem):
        """Generate XML parameters info specific to a peak-calling analysis."""

        ChipseqAnalysis.editParams_XML(self, xmlElem)
        
        paramsNode = ET.SubElement(xmlElem, 'PeakCallingParams')

        # Write options specific to peak calling and profile generation; e.g.
        # profile generation parameters, and MACS parameters...

        # MACS params:
        macsParamsNode = ET.SubElement(paramsNode, 'MACS_params')
        macsParamsNode.text = self.macsOptionsStr

        # Profile generation params:
        profParamsNode = ET.SubElement(paramsNode, 'ProfileParams')
        for profile in self.profileGen.getAllProfs():
            profile.editXMLelem(profParamsNode)

    def writePeakCallsIGV_XML(self, parentNode, baseDataDir2urlBase,
                              nameForXML):
        """Generates IGV data registry XML output for this analysis'
        peak calls."""

        currParamsNode = ET.SubElement(parentNode, "Resource")
        currParamsNode.set("name", nameForXML)

        # Determine BED file URL:
        bedLocAbs = os.path.abspath(self.peaksBedFilename)
        urlLoc = bedLocAbs.replace(baseDataDir2urlBase[0],
                                   baseDataDir2urlBase[1])

        # Record the BED URL location in the XML:
        currParamsNode.set("path", urlLoc)


class SaturationAna(ChipseqAnalysis):
    """Performs peak-calling saturation analysis on a single chip-seq
    experiment."""

    def __init__(self, expt, name, outDir, nReps, fcThreshs,
                 fracsSampled, genome="hg19",
                 macsOptionsStr = "", regStudied="",
                 genomeFile="/home/tom/Work/data/hg19/chromInfo.txt",
                 qValThresh = 5):
        ChipseqAnalysis.__init__(self, expt, name, outDir)

        self.qValThresh = qValThresh

        self.genome = genome

        self.macsOptionsStr = macsOptionsStr

        # Set up dictionaries to allow indexing into the matrix
        # recording the number of peaks called for each replicate,
        # fold-change threshold and fraction sampled...
        self.rep2matrixIdx = {}
        for rep in range(1, nReps + 1):
            self.rep2matrixIdx[rep] = rep-1

        self.fcThresh2matrixIdx = {}
        for fcThreshIdx in range(len(fcThreshs)):
            fcThresh = fcThreshs[fcThreshIdx]
            self.fcThresh2matrixIdx[fcThresh] = fcThreshIdx

        self.fracSampled2matrixIdx = {}
        for fracSampledIdx in range(len(fracsSampled)):
            fracSampled = fracsSampled[fracSampledIdx]
            self.fracSampled2matrixIdx[fracSampled] = fracSampledIdx

        # A 3D matrix recording the number of peaks-called for each
        # replicate, fcThresh, fracSampled combination. All values start as
        # minus one to indicate they have not been set.
        self.nPeaksMatrix = numpy.zeros((len(self.fcThresh2matrixIdx),
                                         len(self.fracSampled2matrixIdx),
                                         len(self.rep2matrixIdx))) - 1

        # The minimum saturation XXX statistic, which provides a single
        # summary statistic indicating the degree of saturation for
        # this experiment.
        self.mser = -1

        # A string specifying the chromosomal coordinate range of the region
        # to conduct peak-calling in. This can be used to do a smaller,
        # faster analysis that nevertheless should produce the same overall
        # results:
        self.regStudied = regStudied

        self.genomeFile = genomeFile

    def getNumPeaksCalled(self, repNum, fcThresh, fracSampled):
        """Retrieves the number of peaks called for the specified, replicate,
        fcCutoff, and fracSampled values."""

        repNumIdx = self.rep2matrixIdx[repNum]
        fcThreshIdx = self.fcThresh2matrixIdx[fcThresh]
        fracSampledIdx = self.fracSampled2matrixIdx[fracSampled]
        return self.nPeaksMatrix[fcThreshIdx][fracSampledIdx][repNumIdx]

    def setNumPeaksCalled(self, repNum, fcThresh, fracSampled, nPeaks):
        """Sets the number of peaks called for the specified, replicate,
        fcCutoff, and fracSampled values."""

        repNumIdx = self.rep2matrixIdx[repNum]
        fcThreshIdx = self.fcThresh2matrixIdx[fcThresh]
        fracSampledIdx = self.fracSampled2matrixIdx[fracSampled]
        self.nPeaksMatrix[fcThreshIdx][fracSampledIdx][repNumIdx] = nPeaks

    def run(self):
        """Runs peak calling, calculates MSER and generates saturation plots."""

        print >> sys.stderr, "Running peak-callings and counting peaks called..."
        self.runPeakCallings()
        print >> sys.stderr, "Done."

        print >> sys.stderr, "Calculating MSER..."
        self.calcMSER()
        print >> sys.stderr, "Done."

        print >> sys.stderr, "Generating plots..."
        self.generatePlots()
        print >> sys.stderr, "Done."

    def getReps(self):
        reps = self.rep2matrixIdx.keys()
        reps.sort()
        return reps

    def getFracsSampled(self):
        fracsSampled = self.fracSampled2matrixIdx.keys()
        fracsSampled.sort()
        return fracsSampled

    def getFcThreshs(self):
        fcThreshs = self.fcThresh2matrixIdx.keys()
        fcThreshs.sort()
        return fcThreshs

    def runPeakCallings(self):
        """Use bamToBed to put all reads (or just those from the region
        being studied) into a temporary bed file, which can be randomly
        sampled from."""

        # Set up the bed files that will be sampled from:
        controlReadsTmpBedFilename = utility.makeTempFilename("TmpBedFileControl")
        conversion.bam2bed(self.expt.getControlAln(),
                           controlReadsTmpBedFilename,
                           regionStr=self.regStudied)

        allReadsTmpBedFilename = utility.makeTempFilename("TmpBedFile")
        conversion.bam2bed(self.expt.getLibSeq(ChipseqExpt.chipName).getBam(),
                           allReadsTmpBedFilename, regionStr=self.regStudied)

        currSampleTmpBedFilename = utility.makeTempFilename("TmpBedFileSample")
        currSampleTmpBamFilename = utility.makeTempFilename("TmpBamFileSample")

        currControlTmpBedFilename = utility.makeTempFilename("TmpBedFileControl")
        currControlTmpBamFilename = utility.makeTempFilename("TmpBamFileControl")

        # For each replicate analysis...
        for rep in self.getReps():
            print >> sys.stderr, "  Replicate", rep, "..."
            # For each fraction of reads sampled...
            for fracSampled in self.getFracsSampled():
                print >> sys.stderr, "  Fraction sampled", fracSampled, "..."

                # Randomly sample that fraction of the reads from the chip
                # and control bed files, and convert to bam:
                currSampleTmpBedFile = open(currSampleTmpBedFilename, 'w')
                utility.sampleLines(allReadsTmpBedFilename,
                                    currSampleTmpBedFile,
                                    fracSampled)
                currSampleTmpBedFile.flush()
                currSampleTmpBedFile.close()
                conversion.bed2bam(currSampleTmpBedFilename,
                                   currSampleTmpBamFilename,
                                   self.genomeFile)

                currControlTmpBedFile = open(currControlTmpBedFilename, 'w')
                utility.sampleLines(controlReadsTmpBedFilename,
                                    currControlTmpBedFile,
                                    fracSampled)
                currControlTmpBedFile.flush()
                currControlTmpBedFile.close()
                conversion.bed2bam(currControlTmpBedFilename,
                                   currControlTmpBamFilename,
                                   self.genomeFile)

                # Run peak calling on the resulting bam file using macs...
                #print >> sys.stderr, "Running peak calling with sampled reads", \
                #    currSampleTmpBamFilename, "..."
                chipLibName = \
                    self.expt.getLibSeq(ChipseqExpt.chipName).getName()
                macsOptionsStr = "--nomodel " + self.macsOptionsStr
                cmdlineProgs.run_MACS(currSampleTmpBamFilename,
                                      currControlTmpBamFilename, self.outDir,
                                      name=chipLibName,
                                      macsOptionsStr = macsOptionsStr,
                                      verbose=True)

                # Check that the peaks.xls file exists and quit, reporting an error,
                # if it does not:
                macsPeaksXlsFilename = self.outDir + "/" + chipLibName + "_peaks.xls"
                try:
                    with open(macsPeaksXlsFilename): pass
                except IOError:
                    print 'Error: MACS failed to produce a peaks.xls file. Quitting.'
                    sys.exit(1)

                # For each fold-change cutoff...
                for fcThresh in self.fcThresh2matrixIdx.keys():
                    #print >> sys.stderr, "  Counting peaks for fcThresh", fcThresh, "..."
                    # Detect the number of peaks called at that threshold:
                    allCurrPeaksTrack = \
                        conversion.macs2bedTrack(open(macsPeaksXlsFilename))
                    peaksPassingThresh = []
                    for peak in allCurrPeaksTrack.getAllBEDs():
                        if ((peak.get_fc() > fcThresh) and
                            (peak.get_pVal() > self.qValThresh)):
                            peaksPassingThresh.append(peak)
                    nPeaksCalled = len(peaksPassingThresh)
                    self.setNumPeaksCalled(rep, fcThresh, fracSampled, nPeaksCalled)
                    print >> sys.stderr, "Rep:", rep, "FcThresh:", fcThresh, \
                        "FracSampled:", fracSampled, "nPeaksCalled:", nPeaksCalled
                    #print >> sys.stderr, "nPeaksMatrix:", self.nPeaksMatrix

        try:
            cmdlineProgs.deleteFiles([allReadsTmpBedFilename,
                                      controlReadsTmpBedFilename,
                                      currSampleTmpBedFilename,
                                      currSampleTmpBamFilename,
                                      currControlTmpBedFilename,
                                      currControlTmpBamFilename])
        except Exception, e:
            print >> sys.stderr, "Failed to delete file(s):"
            print >> sys.stderr, e
            print >> sys.stderr, "Continuing anyway."

    def getNreps(self):
        return len(self.rep2matrixIdx.keys())

    def calcMSER(self):
        """FIXME: Implement later: Compute the MSER statistic on the values in
        self.nRepsMatrix. Shouldn't be too hard."""
        pass

    def getMeanPeaksMatrix(self):
        """Calculates mean and stderr of mean for the number of peaks called
        for each fcThresh and fracSampled value considered. Returns a tuples
        of (meansMatrix, stderrsMatrix)."""

        # Generate a numpy meansMatrix and stderrsMatrix, dimensions determined
        # by number of distinct fcThresh and fracSampled values:
        meansMatrix = numpy.zeros((len(self.fcThresh2matrixIdx),
                                   len(self.fracSampled2matrixIdx))) - 1

        stderrsMatrix = numpy.zeros((len(self.fcThresh2matrixIdx),
                                     len(self.fracSampled2matrixIdx))) - 1

        # For each fcThresh, fracSampled considered...
        for fcThresh in self.getFcThreshs():
            for fracSampled in self.getFracsSampled():
                fcThreshIdx = self.fcThresh2matrixIdx[fcThresh]
                fracSampledIdx = self.fracSampled2matrixIdx[fracSampled]

                # Get every replicate's number of peaks called for that pair,
                # by indexing into self.nPeaksMatrix:
                # FIXME: Will this matrix annotation work?:
                numPeaks_AllReps = \
                    self.nPeaksMatrix[fcThreshIdx][fracSampledIdx]

                # Calculate the mean and standard error of that array of data:
                currMean = numpy.mean(numPeaks_AllReps)
                currStderr = numpy.std(numPeaks_AllReps)/math.sqrt(len(numPeaks_AllReps))

                # Record those in the meansMatrix and stderrsMatrix at
                # fcThreshIdx, fracSampledIdx
                meansMatrix[fcThreshIdx][fracSampledIdx] = currMean
                stderrsMatrix[fcThreshIdx][fracSampledIdx] = currStderr

        return (meansMatrix, stderrsMatrix)

    def generatePlots(self):
        """Generates plot(s) visualising the saturation analysis results."""

        # Generate the 2D matrix of mean number of peaks called for each fracSampled, fcCutoff pair:
        (meanPeaksMatrix, stderrPeaksMatrix) = self.getMeanPeaksMatrix()

        # Generate the output plot file name from self.outDir and self.name:
        self.saturationPlot_filename = self.outDir + "/" + self.name + \
            "_saturationPlot.png"

        # Array for plotting fold-change threshold in colour:
        plotColours = ["black", "blue", "green", "orange", "purple", "grey", "red"]

        maxPeakCount = 0

        # If there are no replicates, then just plot the single replicate:
        if self.getNreps() == 1:
            # Add a separate line plot for each fold-change threshold
            # considered:
            plt.cla()
            plotAxes = plt.axes()
            for fcThresh in self.getFcThreshs():
                # Get the mean number of peaks called for all fractionsSampled,
                # for the current fold-change threshold:
                fcThreshIdx = self.fcThresh2matrixIdx[fcThresh]
                meanPeaksCalled = meanPeaksMatrix[fcThreshIdx]

                # Generate a two-column matrix listing fractionSampled and
                # corresponding mean number of peaks called:
                currDataToPlot = numpy.matrix([self.getFracsSampled(),
                                               meanPeaksCalled.tolist()]).transpose()
                maxPeakCount = max([maxPeakCount, currDataToPlot.max()])

                currColour = plotColours[fcThreshIdx % len(plotColours)]
                utility.addPlotToAxis(currDataToPlot, plotAxes, colour=currColour)

            xlabel = "Fraction of reads sampled"
            ylabel = "Number of peaks called after fold-change filter"
            title = "Saturation analysis: " + self.expt.getName()
            plotFormat="png"

            plotAxes.xaxis.set_label_text(xlabel)
            plotAxes.yaxis.set_label_text(ylabel)
            plotAxes.set_ylim(0, maxPeakCount)
            plotAxes.set_title(title)
            plt.savefig(self.saturationPlot_filename, format=plotFormat)
        else:
            # Not implemented yet:
            print >> sys.stderr, "Multiple replicate saturation analysis " + \
                "not full implemented yet."
            sys.exit(1)
