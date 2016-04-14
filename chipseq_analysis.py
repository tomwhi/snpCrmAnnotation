#!/usr/bin/env python

import commands, math, os, random, shutil, string, sys, pdb
import genomics.formats.bed as bed
import genomics.formats.wibIO as wibIO
import cmdlineProgs, bed1_bed2_dists, getSeqFa, list_funcs, motifModule, motifDists, \
    seq_utility, sequence, seqUtils, utility
import numpy
from optparse import OptionParser


class chipseq_dataset:
    """A single chipseq dataset."""

    def __init__(self, peaksFilename, peaks, reportFilename):
        # The raw input data files for the chipseq dataset:
        self.origPeaksFilename = peaksFilename
        self.reportFilename = reportFilename

        # Generate the bed track corresponding to the peaks now. The input
        # file must be in "peaks bed" format.
        # This class assumes that any unwanted peaks have already been
        # filtered out, and that the input peaks file specified only contains
        # the peaks of interest:
        self.peaks_BED_track = peaks

    def getDatasetName(self):
        return self.getPeaksFilename().replace("/", "_")

    def getPeaksFilename(self):
        return self.origPeaksFilename

    def getPeaksFileTerminus(self):
        return self.origPeaksFilename.split("/")[-1]

    def getReportFilename(self):
        return self.reportFilename

    def getPeaksTrack(self):
        return self.peaks_BED_track

    def writePeakSeqs(self, outFile, flankWidth=250):
        """Introduced on Feb 17th; Write out the peak sequences to specified
        file in fasta format."""

        for peak in self.getPeaksTrack().getAllBEDs():
            # Only use peaks whose sequence has been set successfully:
            if (peak.getPeakSequence(flankWidth, self.getSeqGetter()) != None):
                peak.writeFasta(outFile, flankWidth, self.getSeqGetter())

                seqGetter.getSequence(chrLocQueryStr)


class chipseq_analysis(object):
    """An abstract class, encapsulating a single chipseq analysis. The analysis
    could be on a single dataset or on multiple datasets.

    Keeps track of the output directory for the analysis, as well as the
    name of the analysis."""

    def __init__(self, outDir, analysisName):
        """All analyses know about the output directory to which their output
        is to be written => Store this information in the super constructor,
        here:"""
        # I will treat these variables as "protected"; objects of
        # inheriting classes can access them, but external objects cannot
        # (at least they're not supposed to):
        self.outDir = outDir
        self.analysisName = analysisName
        self.beenRun = False # Indicates whether this analysis has been run yet.
        self.seqGetter = None

    def getOutDir(self):
        return self.outDir

    def getName(self):
        return self.analysisName

    # Hack: September 18th 2013: Refactoring to use a single seqGetter object
    # for all analyses genomic sequence retrievals. Previously, this was
    # done by the bed object.
    def addSeqGetter(self, seqGetter):
        self.seqGetter = seqGetter

    def getSeqGetter(self):
        return self.seqGetter


class multifactor_chipseq_analysis(chipseq_analysis):
    """An abstract class encapsulating a multi-factor chipseq analysis.

    Keeps track of the multiple input datasets for the analysis."""

    def __init__(self, outDir, analysisName, chipseq_datasets):
        super(multifactor_chipseq_analysis,
              self).__init__(outDir, analysisName)

        # chipseq_datasets is supposed to be an array of chipseq_dataset
        # objects. I am not checking this precondition but it should be
        # fullfilled anyway:
        self.chipseq_datasets = chipseq_datasets

    # FIXME: Need to decide how to retrieve the output file for the analysis.
    # One option would be to have a method that infers it from the outDir
    # location - e.g.:
    #self.outFile = outDir + "/QC_Output.html"


class singlefactor_chipseq_analysis(chipseq_analysis):
    """An abstract class encapsulating a single-factor chipseq analysis.

    Keeps track of the input dataset for the analysis."""

    def __init__(self, outDir, analysisName, inputDataset):
        super(singlefactor_chipseq_analysis,
              self).__init__(outDir, analysisName)

        assert (isinstance(inputDataset, chipseq_dataset))
        self.inputDataset = inputDataset

    def getInputDataset(self):
        return self.inputDataset


class chipseq_motif_analysis(singlefactor_chipseq_analysis):
    """An abstract class encapsulating a chipseq motif analysis.

    Keeps track of the input motif library used for the analysis."""

    def __init__(self, outDir, analysisName, inputDataset, library,
                 flankWidth=100, nSeqs=1000):
        super(chipseq_motif_analysis,
              self).__init__(outDir, analysisName, inputDataset)

        # At this stage, the motifLibrary must be the name of a file
        # specifying a list of absolute paths to motif files in MEME format:
        assert (isinstance(library, motifModule.motifLibrary))
        self.library = library # Library of motifs
        self.nSeqs = nSeqs # Number of sequences to sample from input dataset.
        self.flankWidth=flankWidth

    def getMotifLibrary(self):
        return self.library


class singlefactor_chipseq_analyser(object):
    """This abstract class keeps track of a set of analyses to be run/already
    run. The individual analysis objects are kept abstract as chipseq_analysis
    objects. The inheriting concrete classes may assume that the
    chipseq_analysis objects are in fact concrete classes - e.g.
    multifactor_chipseq_analysis or singlefactor_chipseq_analysis classes.
    This is possible because of python's duck typing."""
    def __init__(self):
        # Analyses are added iteratively after the object is created:
        self.analyses = {}

    def getAnalyses(self):
        return self.analyses

    def addAnalysis(self, analysis):
        self.analyses[analysis.getName()] = analysis

    def runAnalyses(self):
        for analysis in self.getAnalyses().values():
            if (not analysis.beenRun):
                print >> sys.stderr, "Running analysis", \
                    analysis.getName(), "..."
                analysis.run()
                print >> sys.stderr, "Completed analysis", analysis.getName(), \
                    "."


class multifactor_chipseq_analyser(object):
    """An instance of this class is used in order to run a set of multifactor
    chipseq analyses."""

    def __init__(self, multifactor_analyses):
        pass # Not implemented yet.


class meme_analysis(singlefactor_chipseq_analysis):
    """A MEME motif analysis of a single input chipseq dataset."""

    def __init__(self, outDir, analysisName, inputDataset, nMEME_seqs,
                 flankWidth=30, nMotifs=5, sampleMethod="highScore",
                 minWidth=5, maxWidth=20, bgFile=None, filterRegsFile=None):
        # A MEME analysis requires that sequences have been set for at least
        # some of the input peaks. This is a precondition of this constructor,
        # although I am currently not enforcing it.

        # Set shared variables of parent class:
        super(meme_analysis, self).__init__(outDir, analysisName, inputDataset)

        if (filterRegsFile != None):
            self.regionFilterTrack = bed.BED_Track(open(filterRegsFile))
        else:
            self.regionFilterTrack = None
        self.nMEME_seqs = nMEME_seqs
        self.flankWidth = flankWidth
        self.memeOutdir = self.outDir + "/MEME_Out_" + self.analysisName
        self.meme_output = None
        self.inFilename = None
        # "highScore"
        self.sample = sampleMethod # Approach towards sampling peaks...
        self.nMotifs=nMotifs
        self.minWidth=minWidth # Minimum motif width to consider
        self.maxWidth=maxWidth # Maximum motif width to consider
        self.bgFile=bgFile

    def get_nMEME_seqs(self):
        return self.nMEME_seqs

    def get_nMotifs(self):
        return self.nMotifs

    def get_memeOutdir(self):
        return self.memeOutdir

    def get_meme_output(self):
        return self.meme_output

    def get_inFilename(self):
        return self.inFilename

    def run(self):
        print >> sys.stderr, "Running MEME analysis..."

        # Select the specified number of strongest peaks from the
        # input dataset, and write their sequences out to a temporary fasta
        # file:
        memeInputSeqsFilename = \
            utility.makeTempFilename(self.outDir + "/" +
                                     self.analysisName + "_MEME_inputSeqs")
        self.inFilename = memeInputSeqsFilename
        memeSeqsFile = open(memeInputSeqsFilename, 'w')
        self.makeSubsetFasta(memeSeqsFile)
        memeSeqsFile.flush()
        memeSeqsFile.close()

        # Run MEME on those input sequences:
        # FIXME: At the moment, I have not implemented any mechanism for
        # passing parameters to MEME. If this is implemented later, it will
        # likely result in some optional parameters being added to the function
        # call here:
        seqUtils.run_MEME(self.get_inFilename(), self.get_memeOutdir(),
                          nMotifs=self.get_nMotifs(), maxWidth=self.maxWidth,
                          minWidth=self.minWidth, bgFile=self.bgFile)

        # Parse the output MEME xml file into an object representing the
        # meme output, including the actual motif data:
        self.meme_output = \
            seqUtils.meme_output(open(self.get_memeOutdir() + "/meme.txt"))

        self.beenRun = True

    # Internal (i.e. "private") methods:
    def makeSubsetFasta(self, outFile):
        """Retrieve the strongest peaks (number specified by attribute of this
        class), and writes their sequences out to a file in fasta format."""

        if (self.regionFilterTrack == None):
            samplePeaks = self.getInputDataset().getPeaksTrack()
        else:
            samplePeaks = self.filterPeaks()

        strongestPeaks = samplePeaks.sampleBEDs(self.nMEME_seqs,
                                                sample=self.sample)

        for peak in strongestPeaks:
            # Only use peaks whose sequence has been set successfully:
            if (peak.getPeakSequence(self.flankWidth, self.getSeqGetter()) != None):
                peak.writeFasta(outFile, self.flankWidth, self.getSeqGetter())

    def filterPeaks(self):
        """Returns a bed track in which each bed item has been filtered to
        make sure it is far (>2kb) away from the closest item in
        self.regionFilterTrack."""

        filteredPeaks = []
        allPeaks = self.getInputDataset().getPeaksTrack().getAllBEDs()
        for peak in allPeaks:
            closestFeature = self.regionFilterTrack.get_closest_item(peak)

            if (closestFeature == None or
                (abs(bed.get_disp(closestFeature, peak)) > 2000)):
                filteredPeaks.append(peak)
        filteredTrack = bed.BED_Track(filteredPeaks)
        return filteredTrack


class xxMotif_analysis(singlefactor_chipseq_analysis):
    """An xxMotif analysis of a single input chipseq dataset."""

    # June 19th 2012: NOTE: I have generated this by copying and pasting
    # and modifying code from meme_analysis. Ideally I should use other
    # methods for code re-use.

    def __init__(self, outDir, analysisName, inputDataset, nXXMotif_seqs,
                 flankWidth=30, sampleMethod="highScore"):
        # Set shared variables of parent class:
        super(xxMotif_analysis, self).__init__(outDir, analysisName, inputDataset)

        self.nXXMotif_seqs = nXXMotif_seqs
        self.flankWidth = flankWidth
        self.xxMotifOutdir = self.outDir + "/XXMotif_Out_" + self.analysisName
        self.xxMotif_output = None
        self.inFilename = None
        self.inputName = None

        self.sample=sampleMethod # Approach towards sampling peaks...

    def get_nXXMotif_seqs(self):
        return self.nXXMotif_seqs

    def get_nMotifs(self):
        return self.nMotifs

    def get_xxMotifOutdir(self):
        return self.xxMotifOutdir

    def get_xxMotif_output(self):
        return self.xxMotif_output

    def get_inFilename(self):
        return self.inFilename

    def run(self):
        print >> sys.stderr, "Running XXMotif analysis..."

        # Select the specified number of strongest peaks from the
        # input dataset, and write their sequences out to a temporary fasta
        # file:
        xxMotifInputSeqsFilename = \
            utility.makeTempFilename(self.outDir + "/" +
                                     self.analysisName + "_XXMotif_inputSeqs")
        self.inFilename = xxMotifInputSeqsFilename
        xxMotifSeqsFile = open(xxMotifInputSeqsFilename, 'w')
        self.makeSubsetFasta(xxMotifSeqsFile)
        xxMotifSeqsFile.flush()
        xxMotifSeqsFile.close()

        # Run XXMotif on those input sequences:
        seqUtils.run_XXMotif(self.get_inFilename(), self.get_xxMotifOutdir())

        inputName=xxMotifInputSeqsFilename.split("/")[-1].split(".")[0]

        self.inputName = inputName

        # Parse the output XXMotif xml file into an object representing the
        # xxMotif output, including the actual motif data:
        self.xxMotif_output = \
            seqUtils.xxMotif_output(self.get_xxMotifOutdir() + "/" + inputName + "_pos_hf.pwm")

        self.beenRun = True

    # Internal (i.e. "private") methods:
    def makeSubsetFasta(self, outFile):
        """Retrieve the strongest peaks (number specified by attribute of this
        class), and writes their sequences out to a file in fasta format."""

        samplePeaks = self.getInputDataset().getPeaksTrack()

        strongestPeaks = samplePeaks.sampleBEDs(self.nXXMotif_seqs,
                                                sample=self.sample)

        for peak in strongestPeaks:
            # Only use peaks whose sequence has been set successfully:
            if (peak.getPeakSequence(self.flankWidth, self.getSeqGetter()) != None):
                peak.writeFasta(outFile, self.flankWidth, self.getSeqGetter())


class spamo_chipseq_analysis(chipseq_motif_analysis):
    """SpaMo analysis of a single input chipseq dataset."""

    def __init__(self, outDir, analysisName, inputDataset, motifLibrary,
                 primaryMotif, flankWidth=250, nSeqs=10000, edgeSpace=150,
                 nIntervals=40, intervalSize=1, eValThresh=0.01, fimoThresh=7,
                 simThresh=0.25, simMetric="Sequence"):
        super(spamo_chipseq_analysis, self).__init__(outDir, analysisName,
                                                     inputDataset,
                                                     motifLibrary,
                                                     flankWidth=flankWidth,
                                                     nSeqs=nSeqs)

#        # Make a parent "spamoResults" folder into which this analysis'
#        # output will be written:
#        self.outDir = outDir + "/spamoResults"
        utility.makeDir(self.getOutDir())

        # The "primary" motif in the context of the SpaMo algorithm:
        self.primaryMotif = primaryMotif

        # SpaMo-specific parameters:
        self.edgeSpace = edgeSpace
        self.intervalSize = intervalSize
        self.nIntervals = nIntervals
        self.eValThresh = eValThresh
        self.fimoThresh = fimoThresh
        self.simThresh = simThresh
        self.simMetric = simMetric

        # Generate a list of tuples representing the intervals to
        # be analysed, given the specified number of intervals and interval
        # size:
        self.intervals = None
        self.generateIntervals()

        # The regions analysed by spamo. Will get set during run():
        self.spamoRegions = None

        # Establish a dictionary "self.results" to store the displacement
        # distribution information for each secondary motif.
        # Each key is the name of a secondary motif that appears in at least
        # one genomic region (empty to start with).
        # Each value is a motifDists.dispResults object showing the distribution
        # of displacements (together with intervals and p-values) for that
        # secondary motif.
        # NOTE: A "back of the envelope" calculation indicates this data
        # structure should typically consume less than 100MB, and will contain
        # all of the spamo results. Presently I deem this acceptable, although
        # the data structure should probably be deleted prior to "pickling".
        self.results = {}

        # Establish a list of dispResults objects that will show the strongest
        # representative results for each "cluster" of results:
        self.nonRedunResults = []

    def getEval(self, motifName):
        """Returns a floating point number equal to the E-value for the
        results associated with the motif of the specified name. E-value
        is calculated by multiplying the result's E-value (after correcting
        for the number of intervals tested) by the number of distinct motifs
        considered."""
        assert(self.results.has_key(motifName))
        nMotifsScanned = self.spamoRegions.get_nMotifsScanned()
        return nMotifsScanned*(self.results[motifName].getEval())

    def getFilteredResults(self):
        """Returns a list of non-redundant dispResults objects that have been
        filtered to exclude those results whose E-value (after correcting for
        the number of motifs tested) is less than the threshold set for this
        object. Assumes that redunReduce() has been run."""
        return filter(lambda result: \
                          (self.getEval(result.secMotif.getName()) <
                           self.eValThresh),
                      self.nonRedunResults)

    def redunReduce(self, simThresh):
        """This method looks at all the results for the various secondary
        motifs, decides which are redundant with which, and stores that
        information in two data structures:
        1) The dispResults objects referenced within self.results
        2) self.nonRedunResults"""

        # Generate an array of all the results (from self.results), sorted by
        # result significance...
        allResults = self.results.values()
        allResults.sort(motifDists.cmpMotifResults)

        # Iterate through all possible result comparisons...
        results1_idx = 0
        while (results1_idx < (len(allResults) - 1)):
            results1 = allResults[results1_idx]

            # Compare the current motifResults against all of the less
            # significant motifResults, and mark those (less significant)
            # results as excluded if they are sufficiently similar to the
            # current one...
            results2_idx = results1_idx + 1
            while (results2_idx < len(allResults)):
                results2 = allResults[results2_idx]

                # See whether the current motifResults and the current less
                # significant motifResults are similar:
                results1results2_sim = results1.calcSim(results2,
                                                        metric=self.simMetric)

                if (results1results2_sim > simThresh):
                    # They are similar => mark the less significant one as
                    # redundant:
                    results1.redundantLessSig.append(results2)
                    results2.declareRedundant()
                results2_idx = results2_idx + 1
            results1_idx = results1_idx + 1

        # Create a list of the non-redundant motif results:
        self.nonRedunResults = filter(lambda mResults: \
                                          not(mResults.isRedundant()), \
                                          allResults)
        self.nonRedunResults.sort(motifDists.cmpMotifResults)

    def run(self):
        """This method does pretty much the same thing as the "main()" method
        of motifDists.py..."""

        # Generate a temporary sequences file of the genomic sequences
        # surrounding the peak summits...
        peakSeqsFilename = \
            utility.makeTempFilename(self.outDir + "/spamo_peakSeqs")
        peakSeqsFile = open(peakSeqsFilename, 'w')

        # Randomly select the specified number of sequences from the input
        # dataset, unless -1 is specified, in which case use all sequences:
        sampledPeaks = self.inputDataset.getPeaksTrack().sampleBEDs(self.nSeqs)

        # Generate a dictionary linking coords to peaks; each key is a string
        # representing the coords of the peak to be scanned, and each value
        # is a peak_BED_line representing the region to be scanned...
        adjustedPeaks = {}
        for peak in sampledPeaks:
            # Convert the peak coordinates so that the start and end of the
            # resulting peaks are each exactly "flankWidth" away from the peak
            # summit:
            adjustedPeak = peak.copy()
            adjustedPeak.setStart(adjustedPeak.get_summit() - self.flankWidth)
            adjustedPeak.setEnd(adjustedPeak.get_summit() + self.flankWidth)
            adjustedPeaks[adjustedPeak.getBedCoordStr()] = adjustedPeak

        # Write the peak sequences out to a temporary fasta file, in
        # preparation for the initial motif scan...
        for peak in adjustedPeaks.values():
            peakSeq = peak.getPeakSequence(self.flankWidth, self.getSeqGetter())
            if (peakSeq != None):
                peak.writeFasta(peakSeqsFile, self.flankWidth, self.getSeqGetter())
        peakSeqsFile.flush()
        peakSeqsFile.close()

        # Generate trimmed, sequences centred on self.primaryMotif...

        # Write the primary motif to a temporary file first:
        tmpPrimaryMotifFilename = utility.makeTempFilename(self.outDir + "/Tmp_MEME_primaryFile_SpaMo", fileSuffix = ".meme")
        tmpPrimaryMotifFile = open(tmpPrimaryMotifFilename, 'w')
        self.primaryMotif.writeToMEME(tmpPrimaryMotifFile)
        tmpPrimaryMotifFile.flush()
        tmpPrimaryMotifFile.close()

        # Generate the centred sequences:
        (centredSeqs, primaryHits) = motifDists.centreSeqs(peakSeqsFilename,
                                                           tmpPrimaryMotifFilename,
                                                           self.edgeSpace,
                                                           self.fimoThresh)

        # Generate a "motifDists.spamoRegs" object, representing the peak
        # sequence regions that have each been scanned by various motifs:
        self.spamoRegions = motifDists.spamoRegs(adjustedPeaks,
                                                 self.primaryMotif,
                                                 primaryHits,
                                                 self.edgeSpace)

        # Write the centred sequences to a second temporary file...
        centredSeqsFilename = \
            utility.makeTempFilename(self.outDir + "/centredSeqs_spamo", \
                                         fileSuffix=".fa")
        centredSeqsFile = open(centredSeqsFilename, 'w')
        seqUtils.writeFasta(centredSeqs, centredSeqsFile)
        centredSeqsFile.flush()
        centredSeqsFile.close()

        # Filter redundant sequences:
        filteredSeqsFilename = \
            utility.makeTempFilename(self.outDir + "/filteredSeqs_spamo", \
                                         fileSuffix=".fa")
        minHammingDist = 150 # Hard-coded; copied from motifDists.main()
        cmd = "redunReduce " + centredSeqsFilename + " " + \
            str(minHammingDist) + " > " + filteredSeqsFilename
        print >> sys.stderr, \
            "Filtering out redundant trimmed, centered sequences; running"\
            + " command " + cmd + "..."
        cmdresult = commands.getstatusoutput(cmd)
        print >> sys.stderr, "Finished:", cmdresult
        filteredSeqs = seqUtils.parseFasta(open(filteredSeqsFilename))

        # For each secondary motif in the library of motifs...
        for secMotifName in self.library.getMotifs():
            print >> sys.stderr, "Analysing spacings for", secMotifName
            secMotif = self.library.getMotif(secMotifName)

            # Write the motif to a temporary file:
            tmpSecMotifFilename = utility.makeTempFilename(self.outDir + "/Tmp_MEME_secFile_SpaMo", fileSuffix = ".meme", check=False)
            tmpSecMotifFile = open(tmpSecMotifFilename, 'w')
            secMotif.writeToMEME(tmpSecMotifFile)
            tmpSecMotifFile.flush()
            tmpSecMotifFile.close()

            print >> sys.stderr, "Analysing secondary motif,", \
                tmpSecMotifFilename, "..."

            # Get the secondary motif hits:
            # NOTE: The "excludedRegion" here was poorly engineered in the
            # original implementation (in motifDists.py). Making it a bit
            # clearer here:
            secondaryHits = seqUtils.getBestHits(filteredSeqsFilename, tmpSecMotifFilename, self.fimoThresh, excludedRegion=[(self.edgeSpace + 1), (self.edgeSpace + 1 + self.primaryMotif.getWidth())])

            # After running a given secondary motif scan, convert the
            # coorindates to be relative to the *original* genomic regions
            # scanned, and then annotate the spamoRegions with the resulting
            # motif hits...
            adjustedSecHits = adjustSecHitLocs(secondaryHits, primaryHits,
                                               self.edgeSpace)
            self.spamoRegions.addSecondaryHits(secMotif, adjustedSecHits)

            # Calculate the distribution of displacements for the current motif:
            self.results[secMotifName] = self.spamoRegions.calcDisps(secMotif)

        # At this point, every secondary motif scan has been performed, and we
        # have a "spamoRegs" object representing all regions scanned by all
        # secondary motifs + the primary motif.

        # Calculate intervals of enrichment for all of the results...
        for result in self.results.values():
            result.calcSigIntervals(self.intervals)

        # Run redundancy reduction over all the motifs:
        self.redunReduce(self.simThresh)

    def writeHtml(self, wibRdrTups=[]):
        """Generates html output for this spamo analysis. Writes all html
        to the specified outFile handle. Generates figures in the specified
        (string) output directory. wibRdrs is an optional list
        of wibIO.wibReader objects, which will be used to output wiggle
        data for each of the results found."""

        outFile = open(self.getOutDir() + "/Output.html", 'w')

        # Write out a header line to the output html file, stating the name of
        # the dataset analysed:
        datasetName = self.inputDataset.getPeaksFilename().split('/')[-1]
        header = '<title>SpaMo analysis for ' + datasetName + '</title>'
        outFile.write(header + '<br/> \n')

        # Generate a sequence logo of the primary motif considered, and
        # generate a link to the resulting file, in the output dir and file
        # respectively...
        primaryMotifLogoPathprefix = self.outDir + "/Spamo_primaryLogo"
        primaryMotifLogoFilename = primaryMotifLogoPathprefix + ".png"
        self.primaryMotif.makeSeqLogo(primaryMotifLogoPathprefix,
                                      format = "png")

        logoAbsPath = os.path.normpath(\
            os.path.abspath(primaryMotifLogoFilename))
        outDirAbsPath = os.path.normpath(os.path.abspath(\
                self.getOutDir()))
        logoRelPath = utility.getRelPath(outDirAbsPath, logoAbsPath)

        logoOutput = '<font size="7">Primary motif:</font><br><img src="' + \
            logoRelPath + '" height="150" alt="MEME motif 1">' + \
            '<br><font size="5">' + self.primaryMotif.getName() + '<font><br>'
        outFile.write(logoOutput + '<br/> \n')

        # Write out the start of an html table:
        tableStart = "<center><big>Results for secondary motifs:</big><br>" + \
            "<table border=\"1\">\n"
        outFile.write(tableStart)
        outFile.flush()

        # Write output for each of the siginficant and non-redundant results...
        for result in self.getFilteredResults():
            # Generate a sequence logo for the secondary motif for this results
            # object in the output directory:
            currSecMotif = result.secMotif
            secMotifLogoPathprefix = \
                self.outDir + "/Spamo_secLogo_" + currSecMotif.getName()
            secMotifLogoFilename = secMotifLogoPathprefix + ".png"
            currSecMotif.makeSeqLogo(secMotifLogoPathprefix,
                                     format = "png")

            # Generate spacing histograms for the current displacement results:
            (ssHistFilename, osHistFilename) = \
                result.makeHists(self.outDir, prefix=currSecMotif.getName())

            # Generate a string representing the statistics for the most
            # significant interval:
            mostSigInterval = result.getMostSigInterval()
            strandStr = "same strand"
            if (not mostSigInterval[2]):
                strandStr = "opp strand"
            sigIntervalStr = str(mostSigInterval[:2]) + "<br>" + strandStr + \
                 "<br>E = %0.3g" % self.getEval(currSecMotif.getName())

            # Start generating the html for the first row for the current
            # results:
            #, and make a link to it in the output file.
            # Make two adjacent links to the hists generated above.
            # Print out the E-value, strand, and bp position of the most
            # significant interval.
            logoAbsPath = os.path.normpath(\
                os.path.abspath(secMotifLogoFilename))
            outDirAbsPath = os.path.normpath(os.path.abspath(\
                    self.getOutDir()))
            logoRelPath = utility.getRelPath(outDirAbsPath, logoAbsPath)

            ssHistAbsPath = os.path.normpath(\
                os.path.abspath(ssHistFilename))
            outDirAbsPath = os.path.normpath(os.path.abspath(\
                    self.getOutDir()))
            ssHistRelPath = utility.getRelPath(outDirAbsPath, ssHistAbsPath)

            osHistAbsPath = os.path.normpath(\
                os.path.abspath(osHistFilename))
            outDirAbsPath = os.path.normpath(os.path.abspath(\
                    self.getOutDir()))
            osHistRelPath = utility.getRelPath(outDirAbsPath, osHistAbsPath)

            currRow = '<tr><td><img src="./' + \
                logoRelPath + \
                '" height="150"><br>' + currSecMotif.getName() + \
                '</img>' + \
                '</td><td><img src="./' + \
                ssHistRelPath + \
                '" height="150"></td><td><img src="./' + \
                osHistRelPath + \
                '" height="150"></td><td>' + \
                sigIntervalStr + \
                '</td>'

            outFile.write(currRow)
            outFile.flush()

            # For retrieving secondary hit genomic sequences:
            seqGetter = getSeqFa.seqRetriever(\
                self.inputDataset.getPeaksTrack().getRefGenomeFilename())

            # Generate wiggle-IC plot output for each of the specified
            # wiggle directories...
            tupIdx = 0
            while tupIdx < len(wibRdrTups):
                tup = wibRdrTups[tupIdx]
                # NOTE: rdrData[0] is a wibReader object, and rdrData[1] is an
                # int specifying the flankWidth to use when drawing profiles.
                # rdrData[2] is a boolean specifying whether to retrieve
                # secondary motif hit locations or primary hit locations.
                rdr = tup[0]
                flankWidth = tup[1]
                getSec = tup[2]

                # Obtain a set of aligned regions corresponding to the
                # positions of the secondary motifs in the sequences exhibiting
                # the most enriched motif spacing:
                sigRegs = result.getSigHits(flankWidth=flankWidth,
                                            getSec=getSec)

                # Obtain the set of sequences corresponding to those
                # regions:
                hitSeqs = []
                for hitReg in sigRegs.getAllBEDs():
                    currSeq = None
                    if (hitReg.getStrand() == "+"):
                        currSeq = seqGetter.getSequence(\
                            hitReg.getBedCoordStr()[3:])
                    else:
                        assert (hitReg.getStrand() == "-")
                        currSeq = seqGetter.getSequence(\
                            hitReg.getBedCoordStr()[3:], revComp=True)
                    hitSeqs.append(currSeq)

                # Obtain a set of aligned regions for the non-enriched regions
                # too:
                nonsigRegs = result.getNonsigHits(flankWidth=flankWidth,
                                                     getSec=getSec)

                # Obtain the set of sequences corresponding to those
                # regions:
                nonsigHitSeqs = []
                for hitReg in nonsigRegs.getAllBEDs():
                    currSeq = None
                    if (hitReg.getStrand() == "+"):
                        currSeq = seqGetter.getSequence(\
                            hitReg.getBedCoordStr()[3:])
                    else:
                        assert (hitReg.getStrand() == "-")
                        currSeq = seqGetter.getSequence(\
                            hitReg.getBedCoordStr()[3:], revComp=True)
                    nonsigHitSeqs.append(currSeq)

                # Obtain a mean wiggle profile array corresponding to
                # the regions spanned by the hits for the regions
                # in the most significant interval:
                sigHitWigProf = aligned_wiggle_analysis.calcAvgWigProfile(\
                    sigRegs, rdr, sigRegs.getAllBEDs()[0].getLen())

                # Obtain a mean wiggle profile for the non-significant regions
                # too:
                sigBEDs = sigRegs.getAllBEDs()[0].getLen()
                nonsigHitWigProf = \
                    aligned_wiggle_analysis.calcAvgWigProfile(\
                    nonsigRegs, rdr, sigBEDs)

                outFile.write('<td>')
                outFile.flush()

                minWigVal = min(sigHitWigProf.getProf() +
                                nonsigHitWigProf.getProf())
                maxWigVal = max(sigHitWigProf.getProf() +
                                nonsigHitWigProf.getProf())
                ylim = (minWigVal, maxWigVal)

                # Generate IC-wiggle html output for the significant and
                # non-significant hit regions:
                outFile.write(rdr.inDir + '<br>')
                outFile.flush()
                currPrefix = "sigSeqs_" + currSecMotif.getName() + str(tupIdx)
                seqUtils.makeSeqAlnWiggleHtml(outFile, self.outDir,
                                              hitSeqs,
                                              sigHitWigProf,
                                              currPrefix, ylim=ylim)
                currPrefix = "nonsigSeqs_" + currSecMotif.getName() + str(tupIdx)
                seqUtils.makeSeqAlnWiggleHtml(outFile, self.outDir,
                                              nonsigHitSeqs,
                                              nonsigHitWigProf,
                                              currPrefix, ylim=ylim)

                outFile.write('</td>')
                outFile.flush()

                outFile.flush()
                tupIdx = tupIdx + 1

            # End the current row in the table:
            outFile.write('</tr>\n')
            outFile.flush()

            # Write out an html table row showing the spaced motif models for
            # the significant intervals...
            interval2model = result.getMotifModels()
            for interval in interval2model.keys():
                model = interval2model[interval]

                # Generate an appropriate unique meme file name:
                currMotifMemeFilename = self.outDir + "/SpacedModel_" + \
                    result.primMotif.name + "_" + result.secMotif.name + \
                    str(interval[0]) + str(interval[2]) + ".meme"

                currMotifMemeFile = open(currMotifMemeFilename, 'w')

                # Output the spaced motif model to that file:
                model.writeToMEME(currMotifMemeFile)
                currMotifMemeFile.flush()
                currMotifMemeFile.close()

                # Infer the relative path to the file:
                modelFileAbsPath = os.path.normpath(\
                    os.path.abspath(currMotifMemeFilename))
                outDirAbsPath = os.path.normpath(os.path.abspath(\
                    self.getOutDir()))
                modelFileRelPath = utility.getRelPath(outDirAbsPath,
                                                      modelFileAbsPath)

                # Write out that relative path as an href link in the html file
                # at this position:
                rowHtml = '<tr><td><a href="' + modelFileRelPath + \
                    '">Model for interval ' + str(interval) + \
                    '</a></td></tr></br>'
                outFile.write(rowHtml)
                outFile.flush()


            # Write an html table row to the output file, showing the secondary
            # motifs whose results are redundant with the current results,
            # and also linking to a bed file showing the locations of the
            # secondary motif hits (strand-specific) for the most significant
            # spacing interval...
            
            # Generate a bed file of the exact "significant" secondary motif
            # locations:
            secHits = result.getSigHits(flankWidth=0, getSec=True)
            secHitsFilename = self.outDir + "/secHits_" + \
                currSecMotif.getName() + ".bed"
            secHitsFile = open(secHitsFilename, 'w')
            secHits.writeToFile(secHitsFile)

            # Start out the html for the row in the table, with a link
            # to the bed file:
            bedAbsPath = os.path.normpath(\
                os.path.abspath(secHitsFilename))
            outDirAbsPath = os.path.normpath(os.path.abspath(\
                    self.getOutDir()))
            bedRelPath = utility.getRelPath(outDirAbsPath, bedAbsPath)
            rowHtml = '<tr><td><a href="' + bedRelPath + \
                '">Significant secondary hits file</a>' + \
                '<br>Redundant Motifs:<br><table>\n'

            # Generate sequence logos and html for each redundant secondary
            # motif...
            for redundResults in result.redundantLessSig:
                redundMotif = redundResults.secMotif

                # Generate current logo:
                currLogoPathprefix = \
                    self.outDir + "/Spamo_secLogo_" + redundMotif.getName()
                currLogoFilename = currLogoPathprefix + ".png"
                redundMotif.makeSeqLogo(currLogoPathprefix,
                                        format = "png")

                # Generate html for current redundant motif:
                logoAbsPath = os.path.normpath(\
                    os.path.abspath(currLogoFilename))
                outDirAbsPath = os.path.normpath(os.path.abspath(\
                        self.getOutDir()))
                logoRelPath = utility.getRelPath(outDirAbsPath, logoAbsPath)

                rowHtml = rowHtml + '<td><img src="' + \
                    logoRelPath + \
                    '" height="50"><br><font size="1">' + \
                    redundMotif.getName() + \
                    '</font></td>'

            # Write the html table row out:
            outFile.write(rowHtml + '</table></td></tr>')
            outFile.flush()

        outFile.flush()

    def generateIntervals(self):
        """Generate a set of intervals to analyse, based on the number of
        intervals and interval size requested."""

        posIntervals = [[(x*self.intervalSize)+1,
                          (x*self.intervalSize+self.intervalSize)]
                         for x in range(int(self.nIntervals/2.0))]
        negIntervals = map(lambda tup: [-tup[1], -tup[0]], posIntervals)
        negIntervals.reverse()
        ssIntervals = map(lambda tup: (tup[0], tup[1], True),
                          negIntervals + posIntervals)
        osIntervals = map(lambda tup: (tup[0], tup[1], False),
                          negIntervals + posIntervals)
        self.intervals = ssIntervals + osIntervals



def adjustSecHitLocs(secHits, primHits, edgeSpace):
    """This function is used to facilitate spamo analysis in the context
    of my new chipseq analysis suite.

    Returns a copy of secHits, in which every motifOccurence
    has been adjusted such that it's start and position have the following
    amount added to it:
    primHit.getStart() - edgeSpace - 1

    This adjustment has the effect of making the resulting hit locations
    relative to the original sequence coordinates, rather than the trimmed
    motif coordinates."""

    adjustedHits = {}
    for seqName in secHits.keys():
        assert(primHits.has_key(seqName))

        secHit = secHits[seqName][0]
        primHit = primHits[seqName][0]
        newStart = secHit.getStart() + primHit.getStart() - edgeSpace - 1
        newEnd = secHit.getEnd() + primHit.getStart() - edgeSpace - 1

        adjustedHit = secHit.copy()
        adjustedHit.setStart(newStart)
        adjustedHit.setEnd(newEnd)

        adjustedHits[seqName] = adjustedHit

    return adjustedHits


class AMA_analysis(chipseq_motif_analysis):
    """An AMA motif analysis of a single input chipseq dataset."""

    def __init__(self, outDir, analysisName, inputDataset, motifLibrary,
                 nTop=5, flankWidth=30, nSeqs=1000):
        super(AMA_analysis, self).__init__(outDir, analysisName, inputDataset,
                                           motifLibrary, flankWidth=flankWidth,
                                           nSeqs=nSeqs)

        # Start out the dictionary of AMA results (motif name keys, pi_one
        # values):
        self.AMAscore_dict = {}
        self.nTop = nTop # Number of top motifs to store.
        self.topMotifNames = None

    def run(self):
        # Obtain a pi_1 value for each motif on the chipseq sequences...

        # Write the peak sequences to a temporary fasta file...
        amaInputSeqsFilename = \
            utility.makeTempFilename(self.outDir + "/AMA_inputSeqs")
        self.inFilename = amaInputSeqsFilename
        amaSeqsFile = open(amaInputSeqsFilename, 'w')

        # Randomly select the specified number of sequences from the input
        # dataset, unless -1 is specified, in which case use all sequences...
        sampledPeaks = self.inputDataset.getPeaksTrack().sampleBEDs(self.nSeqs)

        # Only consider those peaks whose sequence has been set successfully
        # and that doesn't contain "N"s, as "N"s bias the results of an AMA
        # analysis:
        for peak in sampledPeaks:
            peakSeq = peak.getPeakSequence(self.flankWidth, self.getSeqGetter())
            if ((peakSeq != None) and (peakSeq.find("N") == -1)):
                peak.writeFasta(amaSeqsFile, self.flankWidth, self.getSeqGetter())
        amaSeqsFile.flush()
        amaSeqsFile.close()

        # Obtain a background markov model from the input sequences:
        bgFilename = utility.makeTempFilename(self.outDir + "/bgFile")
        seqUtils.run_fastaGetMarkov(amaInputSeqsFilename, bgFilename)

        # Generate a scratch filename for the AMA output:
        amaOutFilename = \
            utility.makeTempFilename(self.outDir + "/AMA_output_scratch")

        # Generate a scratch filename for the p-values:
        pvalsOutFilename = \
            utility.makeTempFilename(self.outDir + "/AMA_pVals_scratch")

        motifLibrary = self.getMotifLibrary().getMotifs()
        for motifName in motifLibrary.keys():
            # Generate a temporary MEME file for the AMA run:
            tmpMemeFilename = \
                utility.makeTempFilename(self.outDir + "/Tmp_MEME_file_AMA",
                                         fileSuffix = ".meme")
            tmpMemeFile = open(tmpMemeFilename, 'w')
            motifLibrary[motifName].writeToMEME(tmpMemeFile)
            tmpMemeFile.flush()
            tmpMemeFile.close()

            # Run AMA on the AMA input sequences for this motif:
            seqUtils.run_AMA(amaInputSeqsFilename, tmpMemeFilename,
                             bgFilename, amaOutFilename)

            cmdlineProgs.deleteFiles([tmpMemeFilename])

            # Parse the results into an ama_results object:
            amaResults = seqUtils.AMA_output(open(amaOutFilename))

            # Output the p-values to a temp scratch file:
            pValsFile = open(pvalsOutFilename, 'w')
            for pVal in amaResults.getPVals():
                pValsFile.write(str(pVal) + "\n")
            pValsFile.flush()
            pValsFile.close()

            currMotif_pi_1 = utility.calc_pi_1(pvalsOutFilename)

            # Store the ama + qValue analysis results for the current motif:
            self.AMAscore_dict[motifName] = currMotif_pi_1

        # Find and store the top motifs...

        # Generate array of motif names and scores to facilitate sorting:
        amaNamesScores = []
        for motifName in self.AMAscore_dict.keys():
            amaNamesScores.append((motifName, self.AMAscore_dict[motifName]))

        # Sort the array on the AMA score values:
        amaNamesScores.sort(utility.cmpTup2, reverse=True)

        self.topMotifNames = map(lambda tup: tup[0], amaNamesScores[:self.nTop])

        # Generate sequence logos for those top motifs:
        for motif in self.topMotifNames:
            currSeqlogo_filePrefix = self.getOutDir() + "/Logo_" + motif
            self.getMotifLibrary().getMotif(motif).makeSeqLogo(\
                currSeqlogo_filePrefix, format="png")

        self.beenRun = True


class bestHit_scanner:
    """A utility class that facilitates easy scanning of chipseq peak sequences
    using a single pwm."""

    def __init__(self, chipseqData, scanPWM, workDir, ownerName, seqGetter,
                 pThresh=0.0001, nSeqs=1000, flankWidth=100, filterN=False):
        assert(isinstance(chipseqData, chipseq_dataset))
        self.chipseqData = chipseqData
        self.scanPWM = scanPWM # The motif.pwm object for scanning
        self.workDir = workDir # A working directory for temp files
        self.pThresh = pThresh
        self.flankWidth = flankWidth

        # A name for an owner of this scanner, which must be unique cf other
        # owners of other bestHit_scanner objects in the program:
        self.ownerName = ownerName

        # Generate an input sequence filename in the workingDir:
        self.inputSeqsFilename = workDir + "/BestHitScanner_InSeqs_" + \
            ownerName + ".fa"

        # Write the peak genomic DNA sequences out to this file...
        seqsFile = open(self.inputSeqsFilename, 'w')

        # Randomly select the specified number of peaks from the input
        # dataset, unless -1 is specified as nSeqs, in which case use all
        # sequences:
        sampledPeaks = None
        if (nSeqs == -1):
            sampledPeaks = self.chipseqData.getPeaksTrack().getAllBEDs()
        else:
            sampledPeaks = self.chipseqData.getPeaksTrack().sampleBEDs(nSeqs)

        # Private variable used for mapping motif occurrences back to the
        # peaks they came from:
        self.str2peak = {}

        # Write out the sequences for those sampled peaks:
        for peak in sampledPeaks:
            # Only use peaks whose sequence can be retrieved successfully,
            # and which do not contain N if that is specified:
            peakSeq = peak.getPeakSequence(self.flankWidth, seqGetter)
            if ((peakSeq != None) and
                ((peakSeq.find("N") == -1) or (not filterN))):
                peak.writeFasta(seqsFile, self.flankWidth, seqGetter)
                self.str2peak[peak.getBedCoordStr()] = peak
        seqsFile.flush()
        seqsFile.close()

        # Generate input and output MAST filenames in the working directory:
        self.motifFilename = workDir + "/BestHitScanner_Motif.mast"
        self.mastOutputFilename = workDir + "/BestHitScanner_MastOutput.txt"

    @staticmethod
    def getAbsHitCoords(bedObject, motifOcc, flankWidth=0):
        """This static method calculates the absolute genomic coordinates
        of the specified motif occurrence, assuming that the sequence scanned
        corresponded to the genomic coordinates specified by the input BED_line
        object. Returns a BED_line object representing the inferred
        coordinates."""

        # Calculate genomic position of motif:
        motifGenStart = bedObject.getChromStart() + motifOcc.getStart() -1 \
            - flankWidth
        motifGenEnd = bedObject.getChromStart() + motifOcc.getEnd() -1 \
            + flankWidth
        motifChrom = bedObject.getChrom()
        motifStrand = motifOcc.getStrand()

        # Return a BED_line object representing the resulting coordinates...
        tmpBedStr = motifChrom + " " + str(motifGenStart) + " " + \
            str(motifGenEnd) + " X 0 " + motifStrand
        return bed.BED_line(tmpBedStr)

    def getBestHits_genomic(self, nHits=1):
        """Returns a dictionary where each key is a peak object and each
        value is a BED_line object representing the chromosomal coordinates
        of the best motif hit.

        Modified on April 19th 2011:
        This method now takes a number nHits=1 as a parameter. This specifies
        the maximum number of best hits the method should return for each peak.
        """

        # Prepare to generate the required dictionary; Retrieve relative
        # motif hit coordinates:
        bestHits = self.getBestHits(nHits=nHits)

        # The dictionary of peaks->genomicHitLocations:
        peak2hitCoords = {}

        # Determine the absolute genomic coordinates of each peak's best
        # hit...
        for currPeak in bestHits.keys():
            # NOTE: This code is getting messy. The underlying problem is that
            # the BED_line class is a little under-engineered and this
            # application is far outside of the class's original intended
            # application. Not sure whether/how to fix this...

            # NOTE: April 19th 2011: This is getting even messier, as I
            # am now modifying the method to allow it to return multiple
            # hits for each peak sequence scanned.

            currPeakTups = [] # Array of (position, motifOccurence) tuples
            for currHit in bestHits[currPeak]:
                # Generate a bed that represents the genomic coordinates of
                # the region that was originally scanned, for the given peak:
                peakSummit = currPeak.get_summit()
                scannedRegChrom = currPeak.getChrom()
                scannedRegStart = peakSummit - self.getSeqFlankWidth()
                scannedRegEnd = peakSummit + self.getSeqFlankWidth()
                scannedRegStr = scannedRegChrom + " " + str(scannedRegStart) + \
                    " " + str(scannedRegEnd) + " X 0 " + currHit.getStrand()
                peakScannedReg = bed.BED_line(scannedRegStr)
                
                currPeakTups.append((bestHit_scanner.getAbsHitCoords(peakScannedReg, currHit), currHit))
                    
            peak2hitCoords[currPeak] = currPeakTups
        return peak2hitCoords

    def getSeqFlankWidth(self):
        return self.flankWidth

    def setPWM(self, newMotif):
        self.scanPWM = newMotif

    def getBestHits(self, nHits=1):
        # Write the pwm out to the MAST file...
        # NOTE: Don't need to do this, since mast accepts meme format now
        # apparently:
        mastFile = open(self.motifFilename, 'w')
        self.scanPWM.writeToMAST(mastFile)
        mastFile.flush()
        mastFile.close()

        # Scan the peak sequences with the motif...
        seqUtils.run_MAST(self.inputSeqsFilename,
                          self.motifFilename,
                          self.mastOutputFilename, pThresh=self.pThresh)

        # Parse the MAST output...
        bestMastHits = seqUtils.getBestMASThits(open(self.mastOutputFilename),
                                                nHits=nHits)

        # Returned dictionary should have peak objects as keys instead of
        # peak sequence names; do this conversion here:
        peak2MastHits = {}
        for peakSeqName in bestMastHits.keys():
            currPeak = self.str2peak[peakSeqName]
            peak2MastHits[currPeak] = bestMastHits[peakSeqName]

        return peak2MastHits


class aligned_wiggle_analysis(singlefactor_chipseq_analysis):
    """This is an analysis in which wiggle data is retrieved for the peak
    regions of the chipseq_dataset, and the average wiggle profile is
    returned. Use this analysis to do nucleosome depletion analysis, for
    example."""

    def __init__(self, outDir, analysisName, inputDataset, wibDir,
                 flankWidth=1000, nPeaks=10000, sampleMethod="random"):
        super(aligned_wiggle_analysis,
              self).__init__(outDir, analysisName, inputDataset)

        # String identifying directory of binary wiggle data:
        self.wibDir = wibDir

        # Width of region to retrieve either side of each peak summit:
        self.flankWidth = flankWidth

        self.nPeaks = nPeaks
        self.sample = sampleMethod

        self.wibRdr = wibIO.wibReader(wibDir)

        # Generate a BED_Track representing the coordinates of the regions
        # surrounding the peak summits...
        currAlignCoords = []
        if (self.nPeaks == -1):
            # Analyse all peaks:
            analysedPeaks = inputDataset.getPeaksTrack().getAllBEDs()
        else:
            analysedPeaks = \
                inputDataset.getPeaksTrack().sampleBEDs(self.nPeaks,
                                                        sample = self.sample)
        for peak in analysedPeaks:
            summitPos = peak.get_summit()
            regStart = summitPos - flankWidth
            regEnd = summitPos + flankWidth
            regChrom = peak.getChrom()
            strand = "+" # Arbitrarily choose plus strand.
            # Only add the BED to the track if it's start and end are valid:
            if (regStart > 0):
                # FIXME: Don't know how to check regEnd; need chromosome
                # sizes to do that. Not sure if this matters or not...
                algnString = regChrom + " " + str(regStart) + " " + \
                    str(regEnd) + " X 0 " + strand
                currAlgnBED = bed.BED_line(algnString)
                currAlignCoords.append(currAlgnBED)
        self.currAlgnBedTrack = bed.BED_Track(currAlignCoords)
        self.currAvgWigProfile = None

    def getCurrAvgProfile(self):
        return self.currAvgWigProfile

    @staticmethod
    def getAllWigProfiles(algnBedTrack, wibRdr, correctLen):
        allWigProfiles = [] # Will be an N*M matrix of wiggle profiles soon...

        # For each alignment position in the current alignment coordinates...
        for pos in algnBedTrack.getAllBEDs():
            # Retrieve the wiggle data from the wib directory for that genomic
            # span. Observe strandedness of the position when deciding whether
            # to "flip" the profile...

            currProfile = wibRdr.getWig(pos, flipOnStrand=True)

            # If a profile was not successfully recovered, it contains
            # null values, or it isn't the correct length, then ignore it:
            if ((currProfile != None) and
                (not currProfile.__contains__(wibRdr.getNull())) and
                (len(currProfile) == correctLen)):
                allWigProfiles.append(currProfile)
        return allWigProfiles

    @staticmethod
    def calcAvgWigProfile(algnBedTrack, wibRdr, correctLen):
        """Calculates the average wiggle profile of a specified track of data,
        given a wib data reader."""

        # Get all of the profiles:
        allWigProfiles = aligned_wiggle_analysis.getAllWigProfiles(algnBedTrack,
                                                                   wibRdr,
                                                                   correctLen)

        if (len(allWigProfiles) == 0 ):
            # Wiggle data could not be retrieved for any regions => return
            # None to indicate this:
            return None

        # Calculate an averaged wiggle profile from the matrix of wiggle data:
        columnSums = reduce(lambda x, y: utility.listSum(x,y), allWigProfiles)
        columnAvgs = []
        idx = 0
        while (idx < len(columnSums)):
            columnAvgs.append(columnSums[idx]/float(len(allWigProfiles)))
            idx = idx + 1
        return wigProfile(columnAvgs)

    def run(self):
        # Simply calculate the average profile for the current alignment
        # coordinates. In the case of this superclass, this means the
        # summit-flanking regions of the input chipseq dataset, as the
        # alignment coordinates are only set once (in the constructor):
        self.currAvgWigProfile = self.calcAvgWigProfile(self.currAlgnBedTrack,
                                                        self.wibRdr,
                                                        ((self.flankWidth*2)
                                                         + 1))
        self.beenRun = True


class detailedHit_analysis(chipseq_motif_analysis):
    """In this analysis, the genomic sequences for the peaks of the specified
    chipseq dataset are scanned (short flanking regions either side of the
    summit), with a specified library of motifs.

    A file of wiggle data is output for each motif. The wiggle track contains
    one block for each peak, with each block showing the LLR contribution of
    each position in the motif hit. This is useful for visualising the
    potential impact of TFBSs at SNP and SNV sites."""

    def __init__(self, outDir, analysisName, inputDataset, motifLibrary,
                 flankWidth=30, pThresh=0.01, nPeaks=100000, nHits=5):
        super(detailedHit_analysis,
              self).__init__(outDir, analysisName, inputDataset, motifLibrary,
                             flankWidth=flankWidth, nSeqs=nPeaks)

        analysisID = \
            analysisName.replace(" ", "_").replace("-", "_").replace("/", "_")

        # Generate a scanner, which will facilitate the motif scans to identify
        # the best hit in each sequence for a given motif:
        self.scanner = bestHit_scanner(self.inputDataset, None, self.outDir,
                                       analysisID, self.getSeqGetter(),
                                       pThresh=pThresh, nSeqs=self.nSeqs,
                                       flankWidth=flankWidth)
        self.wigFilePrefix = outDir + "/motifHits"
        self.nHits = nHits
        self.pThresh = pThresh

    def run(self):
        # For each motif in the input library...
        for motifName in self.library.getMotifs():
            print >> sys.stderr, "Analysing motif", motifName, "..."

            # Obtain the motif of that name:
            currMotif = self.library.getMotif(motifName)
            motifFileprefix = self.outDir + "/" + motifName
            currMotif.makeSeqLogo(motifFileprefix, format="png")

            # Open a new output wiggle track file for the current motif:
            wigFilename = self.wigFilePrefix + "_" + motifName + ".wig"
            wigFile = open(wigFilename, 'w')
            
            # Write out a track header to that file:
            wigFile.write('track type=wiggle_0 name="Motif hits for ' +
                          motifName +
                          '" description="" visibility=2 itemRgb="On"\n')

            # Report that the current motif is being scanned over the chipseq
            # peaks:
            print >> sys.stderr, "Scanning motif " + motifName + " over " + \
                "input dataset..."

            # Select the current motif for scanning:
            self.scanner.setPWM(currMotif)

            # Scan the chipseq dataset being analysed, with the current motif:
            print >> sys.stderr, "Scanning regions..."
            currBestHit_tups = \
                self.scanner.getBestHits_genomic(nHits=self.nHits)

            # Generate a dictionary linking chromosomes to its contained wiggle
            # profiles:
            motifProfs = {}

            # Dictionary of locations where hits have previously been observed,
            # in order to avoid storing the same hit twice:
            prevHits = {}

            # For each of those hits...
            print >> sys.stderr, "Generating wiggle profiles of motif hits..."
            for tups in currBestHit_tups.values():
                # "tups" is an array of (BED_line, motifOccurence) tuples,
                # with each tuple specifying the location and hit statistics of
                # a single motif hit:
                for currTup in tups:
                    # Process the tuple representing the current motif hit
                    # for this ChIP-seq peak...
                    currHit = currTup[0]
                    currMotOcc = currTup[1]
                    chrom = currHit.getChrom()
                    start = currHit.getChromStart()
                    end = currHit.getChromEnd()
                    strand = currHit.getStrand()

                    assert ((end-start+1) == currMotif.getWidth())

                    # Obtain the DNA sequence of this motif hit:
                    hitSequence = None
                    centredPeakLoc = chrom + ":" + str(start) + "-" + str(end)
                    
                    # FIXME: Here, I convert the chromosome string name as
                    # required. This is currently hard-coded, but should be
                    # made more robust and flexible in the future:
                    chrLocQueryStr = centredPeakLoc[3:]
                    seqGetter = getSeqFa.seqRetriever(\
                        self.inputDataset.getPeaksTrack().getRefGenomeFilename())
                    negstrandHit = (strand == "-")
                    hitSequence = seqGetter.getSequence(chrLocQueryStr,
                                                        revComp=negstrandHit)

                    # Get the score contributions of this motif hit, using that
                    # sequence info; wigProf will be an array of motif score
                    # (LLR) contributions.
                    try:
                        wigProf = currMotif.getScoreProf(hitSequence)
                    except Exception, e:
                        pdb.set_trace()
                        x = 1

                    # Add the motif hit profile to the motif occurrence object:
                    currMotOcc.setScoreProfile(wigProf)

                    # Store the motif hit for printing later, but only if
                    # it hasn't already been seen:
                    if (not prevHits.has_key((chrom, start, end))):
                        prevHits[(chrom,start,end)] = True
                        if (motifProfs.has_key(chrom)):
                            motifProfs[chrom].append(currTup)
                        else:
                            motifProfs[chrom] = [currTup]

            # Scanning is complete for this motif. Sort the motif hit
            # information along the chromosomes in preparation for
            # printing:
            print >> sys.stderr, "Sorting motif hit wiggle profiles and " + \
                "filtering out overlaps..."
            for chrom in motifProfs.keys():
                motifProfs[chrom].sort(cmpTups)

                # Filter out overlapping motif hits for the current chromosome;
                # for a given overlapping pair, retain the hit with the lower
                # p-value of the two...
                prevHitIdx=0
                while (prevHitIdx < len(motifProfs[chrom]) - 1):
                    currHitIdx = prevHitIdx+1
                    currHit = motifProfs[chrom][currHitIdx]
                    prevHit = motifProfs[chrom][prevHitIdx]
                    currHitLoc = currHit[0]
                    prevHitLoc = prevHit[0]
                    currHitMotOcc = currHit[1]
                    prevHitMotOcc = prevHit[1]
                    if (prevHitLoc.overlaps(currHitLoc)):
                        print >> sys.stderr, "Overlap between motif hits at", \
                            prevHitLoc.to_string(), "and", \
                            currHitLoc.to_string(), ", filtering out one."
                        # Current and previous hits overlap =>
                        # Remove the hit that has the worst p-value out of the
                        # two:
                        if (prevHitMotOcc.getPval() > currHitMotOcc.getPval()):
                            # Previous hit is worse than the current hit:
                            del motifProfs[chrom][prevHitIdx]
                        else:
                            # Current hit is worse than previous hit:
                            del motifProfs[chrom][currHitIdx]
                    else:
                        # No overlap between current and previous motif hits =>
                        # move to the next position:
                        prevHitIdx = prevHitIdx + 1

            print >> sys.stderr, "Sorted and filtered."

            # Print out the motif occurrences...
            for chrom in motifProfs.keys():
                for tup in motifProfs[chrom]:
                    motifHit = tup[1]
                    # Calculate motifHit's sequence p-value:
                    nScannedPositions = \
                        (self.flankWidth - currMotif.getWidth() + 1) * 4
                    hitPval = motifHit.getPval()
                    seqPval = 1 - \
                        ((1 - hitPval)**nScannedPositions)

                    # Only report match if it passes the motif hit
                    # p-value threshold:
                    if (hitPval < self.pThresh):
                        # Write out a fixedStep header line, stating the
                        # genomic coordinates from currHit, and setting the
                        # colour according to the strand (blue for positive
                        # strand hit, orange for negative strand hit). The
                        # total motif p-value should also be calculated and
                        # appended to the next line after a comment ("#")
                        # character:
                        wigFile.write("fixedStep chrom=" + chrom + " " +
                                      "start=" +
                                      str(tup[0].getChromStart()) +
                                      " step=1\n")
                        wigFile.write("# Hit seq p-val = " + str(seqPval) +
                                      "\n")
                        
                        # Write wiggle data to the output file, representing
                        # the current profile...
                        motifHit.printScoreProfile(wigFile)

            # Flush output to the wiggle file, and close the file:
            wigFile.flush()
            wigFile.close()


def cmpTups(tup1, tup2):
    """Quick hacky comparator function for comparing tuples with both
    bedItem and motifOccurence information in them."""
    hit1Pos = tup1[0].getChromStart()
    hit2Pos = tup2[0].getChromStart()

    if (hit1Pos < hit2Pos):
        return -1
    elif (hit1Pos == hit2Pos):
        return 0
    else:
        assert (hit1Pos > hit2Pos)
        return 1


class bestHit_wiggle_analysis(chipseq_motif_analysis):
    """This is an analysis in which wiggle data is retrieved for the
    regions surrounding the strongest instance of a given motif, for each
    motif in a library. The average wiggle profile is obtained for each motif
    in the library."""

    # Hack: More mess (September 18th 2013) from refactoring seqGetter/bed
    # use: Have to override the base analysis class's addSeqGetter method
    # to also set that info in the scanner:
    def addSeqGetter(self, seqGetter):
        chipseq_analysis.addSeqGetter(self, seqGetter)

        # Generate a scanner, which will facilitate the motif scans to identify
        # the best hit in each sequence for a given motif:
        self.scanner = bestHit_scanner(self.inputDataset, None, self.outDir,
                                       self.analysisID, self.getSeqGetter(), pThresh=self.pThresh,
                                       nSeqs=self.nSeqs,
                                       flankWidth=self.seqScan_flankWidth)

    def __init__(self, outDir, analysisName, inputDataset, wibDir, motifLibrary,
                 flankWidth=1000, nPeaks=10000, sampleMethod="random",
                 cmpFunc="centreVsFlank", cmpInvert=False, pThresh=0.0001):

        # Call chipseq_motif_analysis parent constructor:
        # FIXME: Hard-coded scanning flankWidth, which has a different meaning
        # to the other flankWidth:
        seqScan_flankWidth=100
        super(bestHit_wiggle_analysis,
              self).__init__(outDir, analysisName, inputDataset, motifLibrary,
                             flankWidth=seqScan_flankWidth, nSeqs=nPeaks)

        analysisID = \
            analysisName.replace(" ", "_").replace("-", "_").replace("/", "_")

        # Changed (18-09-2013): Scanner gets set later now.
        self.analysisID = analysisID
        self.scanner = None
        self.pThresh = pThresh
        self.seqScan_flankWidth = seqScan_flankWidth

        # String identifying directory of binary wiggle data:
        self.wibDir = wibDir
        self.sample = sampleMethod
        self.flankWidth = flankWidth
        self.wibRdr = wibIO.wibReader(wibDir)
        self.nPeaks = nPeaks

        # A BED_Track representing the coordinates of the regions surrounding
        # the best hits of the current motif:
        currAlignCoords = None

        # Ranking of motifs in terms of some comparator on the average wiggle
        # profiles:
        self.motifOrdering = None

        # A dictionary containing one average wiggle profile per motif:
        self.motifAvgWiggles = {}

        # Set the comparator function for comparing average wiggle profiles:
        if (cmpFunc == "centreVsFlank"):
            self.profCmp = centreVsFlank_profCmp

        # Set whether to invert the results of the comparator function:
        self.cmpInvert = cmpInvert

        # Select a set of peaks to analyse:
        self.analysedPeaks = None
        if (nPeaks == -1):
            # Analyse all peaks:
            self.analysedPeaks = inputDataset.getPeaksTrack().getAllBEDs()
        else:
            self.analysedPeaks = \
                inputDataset.getPeaksTrack().sampleBEDs(nPeaks,
                                                        sample = self.sample)

    def getFlankWidth(self):
        return self.flankWidth

    def getWibDir(self):
        return self.wibDir

    def setBestHitTrack(self, motif):
        """Helper function, which gets all the best hits for a specified
        motif."""
        # Identify all the best hit genomic coordinates for this motif...
        self.scanner.setPWM(motif)
        currBestHits = self.scanner.getBestHits_genomic(nHits=1)
        
        # Generate a BED_Track of BED_line objects based on those hits
        # but extended to the flankWidth specified...
        bestHitsExtended = []
        # Remember; keys are peak BED_line objects, values are motif hit
        # BED_line representations:
        # Modified on April 19th: getBestHits_genomic() now returns tuples
        # => Extract the location values from that data structure:
        bestHitLocs = map(lambda tup: tup[0][0], currBestHits.values())
        for currHit in bestHitLocs:
            # Generate a new BED_line object based on the original, but
            # with width extended to self.flankWidth...
            chrom = currHit.getChrom()
            start = currHit.getChromStart() - self.flankWidth
            end = currHit.getChromEnd() + self.flankWidth
            strand = currHit.getStrand()
            
            # Only add the BED to the track if it's start and end are valid:
            if (start > 0):
                # FIXME: Don't know how to check regEnd; need chromosome
                # sizes to do that. Not sure if this matters or not...
                tmpBedStr = chrom + " " + str(start) + " " + \
                    str(end) + " X 0 " + strand
                currHitExtended = bed.BED_line(tmpBedStr)
                
                # Add it to the array, from which we will generate the
                # BED_Track:
                bestHitsExtended.append(currHitExtended)
                
        # Generate a BED_Track object from those data:
        self.currAlgnBedTrack = bed.BED_Track(bestHitsExtended)

    def calcAllProfiles(self):
        """Calculates and stores the average wiggle profile for every motif in
        the motif library."""

        # For each motif in the library...
        for currMotif in self.getMotifLibrary().getMotifs().values():
            # Set the alignments track of this object, with the current motif:
            self.setBestHitTrack(currMotif)

            # Get the average wiggle profile for those aligned regions:
            correctProfLen = currMotif.getWidth() + self.flankWidth*2
            currAvgWigProf = aligned_wiggle_analysis.calcAvgWigProfile(self.currAlgnBedTrack, self.wibRdr, correctProfLen)

            # Disregard a motif if no average profile was created:
            if (currAvgWigProf != None):
                currAvgWigProf.setName(currMotif.getName())

                # Store the resulting profile, indexed by the motif:
                self.motifAvgWiggles[currMotif] = currAvgWigProf

    def run(self):
        print >> sys.stderr, "Calculating profiles using width", \
            self.flankWidth, "..."

        # Calculate all of the motifs' average wiggle profiles:
        self.calcAllProfiles()

        # Determine the ranking of the motifs...

        # Create temporary array of the average profiles, to facilitate
        # sorting:
        sortedWigProfiles = self.motifAvgWiggles.values()

        # Do the actual sorting, and obtain a sorted list of motifs:
        sortedWigProfiles.sort(self.profCmp, reverse=self.cmpInvert)
        self.motifOrdering = map(lambda profile: profile.getName(),
                                 sortedWigProfiles)

        self.beenRun = True


def centreVsFlank_profCmp(prof1, prof2):
    """A comparator function that compares two wigProfile objects
    based on the ratio of the central profile score compared with the
    total profile score."""

    if (prof2.getCentralFrac() < prof1.getCentralFrac()):
        return -1
    elif (prof2.getCentralFrac() == prof1.getCentralFrac()):
        return 0
    else:
        assert (prof2.getCentralFrac() > prof1.getCentralFrac())
        return 1


class wigProfile:
    """A simple class representing a profile of wiggle data with an
    associated name."""
    def __init__(self, profile, name=None, fracMiddle=0.33):
        self.profile = profile
        self.name = name
        self.fracMiddle = fracMiddle

        # Calculate ratio of central regions wrt flanking regions...

        # Size of central region:
        profLen = len(profile)
        centralRegSize = int(round(fracMiddle*profLen))

        # Size of flanking regions:
        flankRegSize = int((profLen-centralRegSize)/2.0)

        centralScores = self.profile[flankRegSize:flankRegSize+centralRegSize]
        centralScoreSum = reduce(lambda score1, score2: score1 + score2,
                                 centralScores)

        self.centralScoreSum = centralScoreSum

        self.totalScoreSum = reduce(lambda score1, score2: score1 + score2,
                                    self.profile)

        # Calculate fraction of score that occurs in the central region:
        pseudo = 0.0001 # FIXME: Hack; Not sure how to calibrate pseudo-count
        self.centralFrac = ((self.centralScoreSum + pseudo)/
                            (float(self.totalScoreSum) + pseudo))

    def setName(self, name):
        self.name = name

    def getName(self):
        return self.name

    def getProf(self):
        return self.profile

    def getProfLen(self):
        return len(self.getProf())

    def getCentralFrac(self):
        return self.centralFrac


class motifCentreing_analysis(chipseq_motif_analysis):
    """Motif centreing analysis for a single factor.

     Inputs:
     - A chipseq_dataset object
     - A library of one or more motifs of intereset

     Outputs:
     - Generates histogram data showing distribution of best motif hits
     over peak span, for all motifs in the input library.
"""

    # Hack: More mess (September 18th 2013) from refactoring seqGetter/bed
    # use: Have to override the base analysis class's addSeqGetter method
    # to also set that info in the scanner:
    def addSeqGetter(self, seqGetter):
        chipseq_analysis.addSeqGetter(self, seqGetter)

        # Generate a scanner, which will facilitate the motif scans to identify
        # the best hit in each sequence for a given motif:
        self.scanner = bestHit_scanner(self.inputDataset, None, self.outDir,
                                       self.analysisID, self.getSeqGetter(), pThresh=self.pThresh,
                                       nSeqs=self.nSeqs,
                                       flankWidth=self.flankWidth)

    def __init__(self, outDir, analysisName, inputDataset, motifLibrary,
                 nTop=5, pThresh=0.0001, flankWidth=100, nSeqs=1000,
                 comparator="binom"):
        super(motifCentreing_analysis,
              self).__init__(outDir, analysisName, inputDataset, motifLibrary,
                             flankWidth=flankWidth, nSeqs=nSeqs)

        self.nTop = nTop # Number of top motifs to store.
        self.topMotifNames = None

        # Generate a new bestHit_scanner with the motif as None and the
        # self.input_dataset as the scannedPeaks object. Specify self.outDir
        # as the "workingDir". This will be used to do the MAST scans:
        analysisID = \
            analysisName.replace(" ", "_").replace("-", "_").replace("/", "_")

        # Changed (18-09-2013): Scanner gets set later now.
        self.analysisID = analysisID
        self.scanner = None
        self.pThresh = pThresh

        self.motifHitDistrs = {}

        self.flankWidth=flankWidth

        # Set the comparator function for comparing hit distributions:
        self.hitCmp = None
        if (comparator == "binom"):
            self.hitCmp = binomMotHitCmp
        elif (comparator == "count"):
            self.hitCmp = countMotHitCmp

    def run(self):
        # For each motif in the input library...
        for currMotif in self.getMotifLibrary().getMotifs().values():
            # Prepare the scanner on that motif:
            self.scanner.setPWM(currMotif)

            # Run the motif scan...
            currBestHits = self.scanner.getBestHits(nHits=1)

            # Generate the underlying array of data for a histogram showing the
            # distribution of bestHit locations...

            # Calculate half width of motif and sequence flankWidth, and
            # min and max motif displacement values:
            halfMotW = int(math.floor(currMotif.getWidth()/2.0))
            seqFlankWidth = self.scanner.getSeqFlankWidth()

            # FIXME: Minimum and maximum displacement values actually should
            # consider motif width. However, in order to generate a heatmap
            # over multiple motifs, I will ignore motif width:
            minDisp = -seqFlankWidth #-(seqFlankWidth - halfMotW)
            maxDisp = seqFlankWidth #seqFlankWidth - halfMotW

            # Generate an array of hit displacements based on the dictionary of
            # currBestHits:
            hitPositions = []
            # April 19th 2011: currBestHits now has an array of hits as each
            # value, rather than a single hit => Modifying accordingly using
            # map function:
            for hit in map(lambda arr: arr[0], currBestHits.values()):
                # Calculate centre position of the motif, wrt the start of the
                # scanned sequence:
                hitCentreFromStart = int(math.floor(hit.getStart() + halfMotW))

                # Calculate position wrt the middle of the scanned sequence:
                seqMiddlePos = seqFlankWidth + 1
                hitDispFromMiddle = hitCentreFromStart - seqMiddlePos

                hitPositions.append(hitDispFromMiddle)

            # Generate the histogram data based on the specified array of hit
            # positioning values, and min/max values:

            # 1 bin per 10bp interval:
            nHistBins = int((maxDisp-minDisp+1)/10.0)
            histArr = numpy.histogram(hitPositions, range=[minDisp,maxDisp+1],
                                      bins=nHistBins)[0]

            # Store the motif displacement histogram array for future reference:
            self.motifHitDistrs[currMotif.name] = \
                HitDistr(histArr, minDisp, maxDisp, currMotif.name)

        # Rank the motifs according to the "centredness" of the histograms...

        # Create temporary array of the distributions, to facilitate
        # sorting:
        sortedMotifHit_Distrs = self.motifHitDistrs.values()

        # Do the actual sorting, and obtain sorted list of motifs:
        sortedMotifHit_Distrs.sort(self.hitCmp)
        self.motifOrdering = map(lambda hitDistr: hitDistr.getMotifName(),
                                 sortedMotifHit_Distrs)

        self.beenRun = True
        # Implement later: Instead, sort the motifs ordering based on p-value as
        # calculated using a binomial distribution expected at random (for
        # number of peaks exhibiting a specified motif hit peak displacement
        # value).


class HitDistr:
    """Represents a distribution of motif hits. When designed (Nov 11, 2010),
    this was only intended to be used in the context of a motif centreing
    analysis."""

    def __init__(self, hitDistr, minVal, maxVal, motifName, fracMiddle=0.25):
        self.hitDistr = hitDistr
        self.minVal = minVal
        self.maxVal = maxVal
        self.motifName = motifName

        # Calculate binomial p-value of central hit enrichment, using specified
        # central region...
        nBins = len(hitDistr)
        nCentralBins = int(round(fracMiddle*nBins))
        # Number of bins before the central interval starts:
        nStartBins = int((nBins-nCentralBins)/2.0)

        #print >> sys.stderr, "TRACE:", hitDistr, nBins, nCentralBins, nStartBins, fracMiddle
        centralBins = hitDistr[nStartBins:nStartBins+nCentralBins]
        nCentralHits = int(reduce(lambda count1, count2: count1 + count2,
                                  centralBins))
        self.nCentralHits = nCentralHits
        nTotalHits = int(reduce(lambda count1, count2: count1 + count2,
                                hitDistr))
        self.pVal = utility.binomial(nCentralHits, nTotalHits,
                                     float(nCentralBins)/nBins)

    def getMotifName(self):
        return self.motifName

    def getHitDistr(self):
        return self.hitDistr

    def getNumHits(self):
        """Returns the total number of sequences with a motif hit."""
        return int(reduce(lambda count1, count2: count1 + count2,
                          self.hitDistr))

    def getPval(self):
        return self.pVal

    def getNcentralHits(self):
        return self.nCentralHits


def binomMotHitCmp(motDistr1, motDistr2):
    """Helper function for comparing two motif distributions. Compares based
    on binomial p-value computed for distributions."""

    if (motDistr2.getPval() < motDistr1.getPval()):
        return 1
    elif (motDistr2.getPval() == motDistr1.getPval()):
        return 0
    else:
        assert (motDistr2.getPval() > motDistr1.getPval())
        return -1


def countMotHitCmp(motDistr1, motDistr2):
    """Helper function for comparing two motif distributions. Compares based
    on number of central hits."""

    return motDistr2.getNcentralHits() - motDistr1.getNcentralHits()
