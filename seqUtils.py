#!/usr/bin/env python

# Description:
# This module contains functions for processing sequence data. For example,
# it contains functions for running MAST and parsing the MAST output to produce
# a list of motif occurences.

import commands, math, os, pdb, pylab, random, pandas, sys #, pysam
import cmdlineProgs, motifModule, utility
import numpy as np


class meme_output:
    """This class encapsulates the statistics resulting from a single run
    of MEME. At the moment the class only represents a few stats, although
    I may expand it soon in the future."""

    def __init__(self, memeFile):
        """Parses the specified MEME xml output file, to generate a new object
        representing the stats from that MEME run."""

        self.motifs = {}

        # Parse all the motif data from the memeFile:
        # FIXME: Not sure if below line is robust to input...
        currLine = memeFile.readline()
        while (currLine != ""):
            try:
                currMotif = motifModule.pwm(memeFile)
                self.motifs[currMotif.getName()] = currMotif
            except ValueError, e:
                # FIXME: DODGY! This is expected if no motif is found, which
                # should occur when the meme file has been exhausted of motifs.
                # However, this flow control is poor...
                pass

            # FIXME: This is hacky: The problem is that the parsing is being
            # done by disparate objects that don't share the same state in
            # terms of current line of the file.
            currLine = memeFile.readline()

    def getMotifs(self):
        return self.motifs

    def getMotifsArr(self):
        """Returns an array of the motifs, sorted from first to last in order
        of discovery."""
        motifNum = 1
        motifsSorted = []
        while (motifNum <= len(self.motifs.values())):
            currKey = str(motifNum)
            currMotif = self.motifs[currKey]
            motifsSorted.append(currMotif)
            motifNum = motifNum + 1
        return motifsSorted


class xxMotif_output:
    """This class encapsulates the statistics resulting from a single run
    of XXmotif. At the moment the class only represents a few stats, although
    I may expand it soon in the future."""

    def __init__(self, xxmotifFilename):
        """Parses the specified XXmotif xml output file, to generate a new
        object representing the stats from that XXmotif run."""

        self.motifs = {}

        try:
            xxmotifFile = open(xxmotifFilename)
        except IOError, e:
            print >> sys.stderr, "WARNING: Could not open file:", xxmotifFilename
            return

        # Parse all the motif data from the xxmotifFile:
        # FIXME: Not sure if below line is robust to input...
        currLine = None
        while (currLine != ""):
            try:
                currMotif = motifModule.pwm(xxmotifFile, dataType="xxMotifFile")
                self.motifs[currMotif.getName()] = currMotif
            except ValueError, e:
                # FIXME: DODGY! This is expected if no motif is found, which
                # should occur when the xxmotif file has been exhausted of motifs.
                # However, this flow control is poor...
                pass

            # FIXME: This is hacky: The problem is that the parsing is being
            # done by disparate objects that don't share the same state in
            # terms of current line of the file.
            currLine = xxmotifFile.readline()

    def getMotifs(self):
        return self.motifs

    def getMotifsArr(self):
        """Returns an array of the motifs, sorted from first to last in order
        of discovery."""
        motifNum = 1
        motifsSorted = []
        while (motifNum <= len(self.motifs.values())):
            currKey = str(motifNum)
            currMotif = self.motifs[currKey]
            motifsSorted.append(currMotif)
            motifNum = motifNum + 1
        return motifsSorted


class AMA_output:
    """This class encapsulates the statistics resulting from a single run
    of AMA"""

    def __init__(self, amaFile):
        """Parses the specified cisml output file, to generate a new object
        representing the stats from that AMA run."""

        self.seqName2pVal = {}

        # NOTE: At the moment I just parse the sequence names and p-values

        # Parse the sequence names and p-values...
        for line in amaFile.xreadlines():
            elems = line.split()
            if ((len(elems) >= 5) and (elems[4][:6] == "pvalue")):
                name = elems[2].split("=")[1].strip('"')
                pVal = float(elems[4].split("=")[1].strip('"'))
                #self.seqStats[name] = pVal
                # IGNORING NAME NOW - Just store the pvalue:
                self.seqName2pVal[name] = pVal

    def getPVals(self):
        return self.seqName2pVal.values()

    def getPeak2pVal(self):
        return self.seqName2pVal


def getBestHits(seqsFilename, motifFilename, fimoThresh, space=0, \
                excludedRegion=None):
    """This function takes a fasta file of sequences and a motif filename as input, and returns a hashtable (with sequence name keys and motifOccurrence values) indicating the best hit of the given motif in each sequence."""

    # FIXME: Use FIMO by default but perhaps in the future introduce an option
    # to use MAST. If I do that, then it might be wise to use an object-oriented
    # approach, with a "scanner" abstract class and "mastScanner" and
    # "FIMO_Scanner" concrete classes (or something along those lines).

    # Create a temporary file to hold the *processed* (i.e. one hit per
    # sequence) FIMO output:
    processedHitsFilename = utility.makeTempFilename("ProcessedHits")
    processedHitsFile = open(processedHitsFilename, 'w')

    # Scan the sequences with the motif, printing the best hits to the
    # temporary file:
    printBestHits(seqsFilename, motifFilename, processedHitsFile, fimoThresh, space=space, excludedRegion=excludedRegion)
    processedHitsFile.flush()
    processedHitsFile.close()

    # Parse the resulting hits, to obtain the hashtable of hits:
    bestHits = parseFimo(open(processedHitsFilename))

    cmdlineProgs.deleteFiles([processedHitsFilename])
    return bestHits


def processFimoOutput(fimoFile, outfile, seqs, space=0, excludedRegion=None):
    """Processes the specified file containing output from fimo, outputting
    the first, best pwm hit location for each sequence that was scanned by
    fimo."""
    currSeqName = None
    currSeqScores = None
    for line in fimoFile:
        elems = line.strip().split("\t")
        if (elems[0][0] != "#"):
            # The line is not a comment line, => Parse it.
            if(len(elems) != 8):
                raise Exception("Invalid line: " + line)
            
        try:
            currHit = parseFimoLine(line)

            if (currHit.getSeqName() == currSeqName):
                # Line represents data for the same sequence as was considered
                # in the previous line. => Add it to the set of scored nucs for
                # this sequence:
                currSeqScores.append(currHit)
            else:
                # Have encountered the start of a new sequence. => Print out
                # information for the previous sequence (if there was any info),
                # and start to consider this new sequence...
                if (currSeqName != None):
                    processSequenceBlock(currSeqName, currSeqScores, outfile, \
                                         seqs, space=space, \
                                         excludedRegion=excludedRegion)

                # Start considering the next sequence...
                currSeqName = currHit.getSeqName()
                currSeqScores = [currHit]

        except ValueError, e:
            print >> sys.stderr, "Invalid fimo output line:\n", line
            raise(e)
            
    # Process the final sequence block, if there was one:
    if (currSeqName != None):
        processSequenceBlock(currSeqName, currSeqScores, outfile, seqs, \
                             space=space, excludedRegion=excludedRegion)
      

def processSequenceBlock(seqName, seqScores, outfile, seqs, space=0, \
                         excludedRegion=None):
    """Finds the single best-scoring motif occurrence in the specified list
    of scores (that must all belong to the same specified sequence name), within
    the specified range of the sequence. Prints that motif occurence, if one
    was found."""

    # PRECONDITION:
    assert(len(seqs[seqName]) > space*2)

    coords = seqName.strip().split(".")[-1]
    chrom = coords.split(":")[0]
    startStop = coords.split(":")[1]
    try:
        seqStart = int(startStop.split("-")[0])
        seqStop = int(startStop.split("-")[1])
    except utility.FormatError, e:
        print >> sys.stderr, "Invalid chromosomal coordinate string:", \
            coords
        sys.exit(1)

    # XXX CONTINUE HERE: Apparent problem: seqScores is too short apparently;
    # There should be a score for each nucleotide position in the sequence!
    # It's because psp_fimo is not printing info for each nucleotide, as I had
    # assumed it would.

    # Filter all hits to *exclude* those that occur closer than <edgeSpace>
    # to the edge of the given sequence, or which overlap with the excluded
    # region:
    # NOTE: Excluded region coordinates are specified as the same for every
    # sequence. Thus this parameter is probably only useful when the sequences
    # are aligned from their start (as is the case when they are aligned on a
    # primary motif instance and trimmed).
    currSeq = seqs[seqName]
    centralHits = filter(lambda hit: (hit.getStart() > space) and \
                         (hit.getEnd() < (len(currSeq) - space)), seqScores)

    # "filteredHits" is the set of central hits that have had the specified
    # region excluded:
    filteredHits = filter((lambda hit: (excludedRegion == None) or \
                          (hit.getEnd() < excludedRegion[0]) or \
                          (hit.getStart() > excludedRegion[1])), centralHits)

    # Find the best hit amongst the remaining hits (resolving ties by
    # randomisation)...

    # Get all hits with the equal maximum score in the current sequence:
    INFINITY = 1000000000
    maxScore = -INFINITY
    for hit in filteredHits:
        if (hit.getScore() > maxScore):
            maxScore = hit.getScore()
    
    bestHits = filter(lambda x: x.getScore() == maxScore, filteredHits)

    # Declare the first of the best hits as the gold standard nucleotide
    # for the current region: FIXME: THAT COMMENT NO LONGER APPLIES;
    # BUT I NEED TO GET makeGoldStandard.py working again before I should
    # delete the comment.

    if (len(bestHits) > 0):
        # A hit (and therefore a best hit) was identified for this sequence
        # block => print it out...

        bestHit = random.choice(bestHits) #bestHits[0]
        # FIXME: The next line of code was used when this function did
        # makeGoldStandard stuff; need to resolve that change in order to
        # get that program to work again.
#    siteStart = seqStart + curr_gs_nuc.getStart()
#    print >> outfile, chrom, siteStart, curr_gs_nuc.getStrand(), \
#        curr_gs_nuc.getScore(), curr_gs_nuc.getMatch(), len(bestHits)
        
        # Print the best hit in FIMO output format:
        print >> outfile, bestHit.toFIMO_String()


def printBestHits(seqsFilename, motifFilename, outFile, fimoThresh, space=0, excludedRegion=None):
    """This function scans the given sequences with the given motif, and prints the output to the specified output file."""

    # Create a temporary filename, for the file to hold the fimo output:
    fimo_output_tmp_filename = utility.makeTempFilename("fimo_output")

    # Run fimo, placing output into a file of that name:
    runPspFimo(motifFilename, seqsFilename, fimo_output_tmp_filename, fimoThresh)

    # Parse the sequence file, so that required sequence information (such as
    # sequence length) is available:
    seqs = parseFasta(open(seqsFilename))

    # Process the resulting motif hits, printing the single best hit in each
    # sequence:
    processFimoOutput(open(fimo_output_tmp_filename), outFile, seqs, \
                      space=space, excludedRegion=excludedRegion)

    # Delete the temporary files:
    cmdlineProgs.deleteFiles([fimo_output_tmp_filename])


def run_tomtom(motifFilename, motifLibrary, outDir):
    """Runs tomtom (assumes that the command 'tomtom' is in the path), and
    specifies the directory name <outDir> as the output location."""

    # FIXME: Need to introduce some way of checking whether the command is
    # successful, and reporting an error otherwise.

    cmd = "tomtom -query " + motifFilename + " -target " + motifLibrary + \
        " -oc " + outDir

    print >> sys.stderr, "Running meme, with command:\n", cmd
    cmdResult = commands.getstatusoutput(cmd)
    print >> sys.stderr, "Command result:", cmdResult


def run_MEME(seqsFilename, outDir, nMotifs=5, minWidth=5, maxWidth=20,
             bgFile=None):
    """Runs MEME (assumes that the command 'meme' is in the path), and specifies
    the directory name <outDir> as the output location."""

    # FIXME: Need to introduce some way of checking whether the command is
    # successful, and reporting an error otherwise.

    cmd = "meme -maxsize 2000000 -nmotifs " + str(nMotifs) + \
        " -revcomp -dna " + seqsFilename + " -o " + outDir + " -minw " + \
        str(minWidth) + " -maxw " + str(maxWidth)

    if (bgFile != None):
        cmd = cmd + " -bfile " + bgFile

    print >> sys.stderr, "Running meme, with command:\n", cmd
    cmdResult = commands.getstatusoutput(cmd)
    print >> sys.stderr, "Command result:", cmdResult


def run_XXMotif(seqsFilename, outDir):
    """Runs XXMotif (assumes that the command 'XXmotif' is in the path), and
    specifies the directory name <outDir> as the output location."""

    cmd = "/home/tom/Work/externalProgs/XXmotif-64Bit-static/XXmotif " + outDir + " " + seqsFilename + " --mops --revcomp --XXmasker"

    print >> sys.stderr, "Running xxmotif, with command:\n", cmd
    cmdResult = commands.getstatusoutput(cmd)
    print >> sys.stderr, "Command result:", cmdResult


def run_fastaGetMarkov(seqsFilename, bgFilename, order=3):
    """Runs fasta-get-markov (assumes that the command is in the path), and
    specifies the directory file <bgFilename> as the output location."""

    cmd = "fasta-get-markov -m " + str(order) + " < " + seqsFilename + " > " + bgFilename

    print >> sys.stderr, "Running fasta-get-markov, with command:\n", cmd
    cmdResult = commands.getstatusoutput(cmd)
    print >> sys.stderr, "Command result:", cmdResult


def run_ceqlogo(matrixFilename, outFilePrefix, format="eps"):
    """Runs fasta-get-markov (assumes that the command is in the path), and
    specifies the directory file <bgFilename> as the output location."""

    cmd = "ceqlogo -t \"\" -d \"\" -x \"\" -i " + matrixFilename + " -m 1 > " + outFilePrefix + ".eps"

    print >> sys.stderr, "Running ceqlogo, with command:\n", cmd
    cmdResult = commands.getstatusoutput(cmd)
    print >> sys.stderr, "Command result:", cmdResult

    if (format == "eps"):
        # No image conversion needed:
        pass

    if (format == "png"):
        # Convert to png...

        cmd = "convert " + outFilePrefix + ".eps " + outFilePrefix + ".png"

        print >> sys.stderr, "Running convert, with command:\n", cmd
        cmdResult = commands.getstatusoutput(cmd)
        print >> sys.stderr, "Command result:", cmdResult

        cmdlineProgs.deleteFiles([outFilePrefix + ".eps "])


def run_AMA(seqsFilename, motifFilename, bgFilename, outputFilename):
    """Runs AMA (assumes that the commands 'ama' is in the path), and
    specifies the file <outputFilename> as the output location."""

    # FIXME: Need to introduce some way of checking whether the command is
    # successful, and reporting an error otherwise.

    cmd = "ama --pvalues " + motifFilename + " " + seqsFilename + " " + \
        bgFilename + " > " + outputFilename

    print >> sys.stderr, "Running AMA, with command:\n", cmd
    cmdResult = commands.getstatusoutput(cmd)
    print >> sys.stderr, "Command result:", cmdResult


def run_MAST(sequenceFilename, motifFilename, outputFilename, pThresh = 0.01):
    """Runs MAST (assumes that the command 'mast' is in the path), and specifies
    <outputFilename> as the output location."""

    cmd = "mast " + motifFilename + " " + sequenceFilename + \
        " -nohtml -hit_list -mt " + str(pThresh) + " > " + outputFilename

    print >> sys.stderr, "Running genome scan, with command:\n", cmd
    cmdResult = commands.getstatusoutput(cmd)
    print >> sys.stderr, "Command result:", cmdResult


def run_strcount(sequenceFilename, outFilename, options = ""):
    """Runs strcount (assumes that the command 'strcount' is in the path), and
    specifies <outFilename> as the output location."""

    cmd = "strcount " + options + " -s " + sequenceFilename + \
        " -o " + outFilename

    print >> sys.stderr, "Running strcount, with command:\n", cmd
    cmdResult = commands.getstatusoutput(cmd)
    print >> sys.stderr, "Command result:", cmdResult


def parseFimoLine(line, fimoType="pspFimo"):
    """Parses a single line of FIMO output, producing a motifOccurence object"""

    elems = line.strip().split()

    # April 15th 2013: Detect if regular fimo or psp_fimo output format,
    # and parse accordingly:
    if fimoType == "pspFimo":
        motifName = elems[0]
        seqName = elems[1]
        start = int(elems[2])
        end = int(elems[3])
        score = float(elems[4])
        pVal = elems[5]
        qVal = elems[6]
        matchedSeq = elems[7]
    
        strand = "+"
        motifStart = start
        motifEnd = end
        if (start > end):
            strand = "-"
            # I define start of the binding site as the lowest genomic
            # coordinate of it:
            motifStart = end
            motifEnd = start
    else:
        assert (fimoType == "fimo")

        motifName = elems[0]
        seqName = elems[1]
        motifStart = int(elems[2])
        motifEnd = int(elems[3])
        strand = elems[4]
        score = float(elems[5])
        pVal = float(elems[6])
        qVal = None
        matchedSeq = elems[7]

    return motifOccurence(seqName, strand, motifStart, motifEnd, score=score, \
                          motifName=motifName, pVal=pVal, qVal=qVal, \
                          matchedSeq=matchedSeq)


# Introduced on June 19th 2013:
def parseFimoBlock(currLine, fimoOutput, fimoType="pspFimo"):
    """Assumes the input fimo file is sorted by motif. Parses all fimo output
    for a single motif. Returns a sequence2pwmHit dictionary.

    "currLine" is the current line of input from the specified open file."""
    seq2pwmHit = {}
    # Indicates when new motif block is encountered in the fimo output file:
    newMotif = False
    while (currLine != "" and newMotif == False):
        elems = currLine.strip().split()
        currMotif = elems[0]
        if (elems[0][0] != "#" and elems[0] != "Motif"):
            # Line is not a comment and hence should be parsed...
            if (len(elems) != 8 and fimoType == "pspFimo"):
                print >> sys.stderr, "Invalid fimo line:", line, \
                    "; wrong number of elems:", len(elems)
                raise Exception()
            try:
                # Create a new motif occurence representing the current line
                # (i.e. current hit):
                currHit = parseFimoLine(currLine, fimoType=fimoType)
                if (seq2pwmHit.has_key(currHit.getSeqName())):
                    seq2pwmHit[currHit.getSeqName()].append(currHit)
                else:
                    seq2pwmHit[currHit.getSeqName()] = [currHit]
                    
            except ValueError, e:
                print >> sys.stderr, "Invalid fimo output line:\n", currLine
                raise(e)

        currLine = fimoOutput.readline()
        if (currLine != "" and currLine.strip().split()[0] != currMotif):
            # The next line has a different motif line => It belongs to a new
            # block:
            newMotif = True

    print >> sys.stderr, \
        "Processed fimo output file block for motif %s." % currMotif

    return (currLine, seq2pwmHit)


def parseFimo(fimoOutput, fimoType="pspFimo"):
    """Parses a specified stream of fimo (text-mode) output, to produce a hash
    data structure where each fasta sequence (in the fimo output) has a
    corresponding entry in the hashtable. Each hashtable entry has the name of
    the sequence as the key, and an array representing the pwm hits (in that
    sequence) as the value."""
    # FIXME: Modified output data structure to be an array of such dictionaries,
    # one per motif block:
    allFimoData = []

    # FIXME: Haven't re-tested this after modifying on June 19th 2013.
    currLine = fimoOutput.readline()
    currBlockData = parseFimoBlock(currLine, fimoOutput, fimoType)
    allFimoData.append(currBlockData)
    while currLine != "":
        (currLine, currBlockData) = parseFimoBlock(currLine, fimoOutput, fimoType)
        allFimoData.append(currBlockData)

    # FIXME: This is broken now too:
    for key in outputHash.keys():
        outputHash[key].sort()

    return outputHash


def runPspFimo(pwmFileName, sequenceFileName, outputFileName, scoreThresh, options=""):
    cmd = "psp_fimo " + options + " --output-score_thresh " + str(scoreThresh) + " " + pwmFileName + " " + sequenceFileName + " > " + outputFileName
    print >> sys.stderr, "Running FIMO, with this command:\n", cmd
    cmdresult = commands.getstatusoutput(cmd)
    print >> sys.stderr, cmdresult


def runFimo(pwmFileName, sequenceFileName, outputFileName, options="", verbose=True):
    cmd = "fimo --text " + options + " " + pwmFileName + " " + sequenceFileName + " > " + outputFileName
    if verbose:
        print >> sys.stderr, "Running FIMO, with this command:\n", cmd
    cmdresult = commands.getstatusoutput(cmd)
    if verbose:
        print >> sys.stderr, cmdresult


def run_GLAM2SCAN(sequenceFilename, motifFilename, outputFilename, nHits = 25):
    """Runs GLAM2SCAN (assumes that the command 'glam2scan' is in the path),
    and specifies <outputFilename> as the output location."""

    cmd = "glam2scan -o " + outputFilename + " -2 -n " + str(nHits) + " n " + \
        motifFilename + " " + sequenceFilename

    #print >> sys.stderr, "Running GLAM2SCAN, with command:\n", cmd
    cmdResult = commands.getstatusoutput(cmd)
    #print >> sys.stderr, "Command result:", cmdResult


class motif:
    """A motif model of TFBS binding. In introduced this on 4th of October 2010,
    to help me with parsing MEME output files. At the moment it just contains
    a few stats, although I could potentially expand this class in the future.
    If I decide to expand the class, I should first check whether there is
    another motif rep class somewhere in my code base, or whether there is
    some such useful code implemented already by a third party."""

    def __init__(self, eVal):
        self.eValue = eVal

    def getEvalue(self):
        return self.eValue


class motifOccurence:
    """A motif occurrence."""
    # FIXME: At some point perhaps I should re-think how I represent motif
    # occurences. An object oriented approach sounds suitable here, as there
    # are "mast" hits and "fimo" hits, which could each implement an
    # abstract class.

    def __init__(self, seqName, strand, start, end, pVal = None,
                 motifName = None, score = None, qVal = None,
                 matchedSeq = None):
        self.seqName = seqName
        if (strand == "+1" or strand == "+"):
            self.strand = "+"
        elif (strand == "-1" or strand == "-"):
            self.strand = "-"
        else:
            raise ValueError("Invalid strand value: \"" + strand + "\"")

        # NOTE: "Start" and "End" indeces are *inclusive* of the motif, and
        # count from 1 as the first base in the sequence (not zero).
        self.start = start
        self.end = end
        self.pVal = pVal
        self.qVal = qVal
        self.score = score
        self.motifName = motifName
        self.score = score
        self.matchedSeq = matchedSeq
        # Introduced on 5th November 2010; A motif.pwm object that this
        # is an occurrence of:
        self.motif = None
        # Introduced on May 27th 2011: The score profile of this motif
        # occurrence:
        self.scoreProfile = None

    def copy(self):
        return motifOccurence(self.seqName, self.strand, self.start, self.end, pVal=self.pVal, motifName=self.motifName, score=self.score, qVal=self.qVal, matchedSeq=self.matchedSeq)

    def getScoreProfile(self):
        return self.scoreProfile

    def setScoreProfile(self, profile):
        self.scoreProfile = profile

    def printScoreProfile(self, outFile):
        if (self.scoreProfile == None):
            print >> outFile, None
        else:
            prof = list(self.scoreProfile)
            if (self.getStrand() == "-"):
                prof.reverse()
            for score in prof:
                outFile.write(str(score) + "\n")

    def setMotif(self, newMotif):
        self.motif = newMotif

    def getMotif(self):
        return self.motif

    def getPval(self):
        return self.pVal

    def getScore(self):
        return self.score

    def getSeqName(self):
        return self.seqName

    def getStart(self):
        return self.start

    def getEnd(self):
        return self.end

    def setStart(self, start):
        self.start = start

    def setEnd(self, end):
        self.end = end

    def getStrand(self):
        return self.strand

    def toFIMO_String(self):
        """Returns a string representing the motif occurrence in the FIMO text
        output format."""
        # FIMO doesn't represent strand explicitly, but does it implicitly
        # by the order of the start and end:
        if (self.getStrand() == "+"):
            motifStart = self.start
            motifEnd = self.end
        else:
            motifStart = self.end
            motifEnd = self.start
            
        return self.motifName + " " + self.seqName + " " + str(motifStart) + \
            " " + str(motifEnd) + " " + str(self.score) + \
            " " + str(self.pVal) + " " + str(self.qVal) + " " + self.matchedSeq

    def toString(self):
        #print >> sys.stderr, str(self.start), str(self.end), str(self.strand),\
        #    str(self.pval)
        return str(self.start) + " " + str(self.end) + " " + str(self.strand) \
            + " " + str(self.pVal) + " " + str(self.motifName) + " " + \
            str(self.score)

    def getMotifName(self):
        return self.motifName

    def getMatchedSeq(self):
        return self.matchedSeq

    def getCentre(self):
        return (self.getStart() + self.getEnd())/2.0

    def sameStrand(self, mot2):
        return self.getStrand() == mot2.getStrand()

    def getWidth(self):
        return self.getEnd() - self.getStart()

    def getDisp(self, mot2):
        """Returns the displacement between this motif occurrence and the
        specified motif occurrence, relative to the strand that this motif
        occurs on. Returns None if the motifs overlap."""

        posStrandDisp = None
        if (mot2.getEnd() < self.getStart()):
            # The end of the second occurrence is to the left of the start
            # of this one => negative displacement:
            posStrandDisp = -(self.getStart() - mot2.getEnd())
        else:
            if (mot2.getStart() <= self.getEnd()):
                # The motifs overlap, => return None:
                return None
            else:
                # The start of the second occurrence is to the right of the
                # end of this one => positive displacement:
                posStrandDisp = (mot2.getStart() - self.getEnd())

        if (self.getStrand() == "+"):
            return posStrandDisp
        else:
            # This motif is on the negative strand, so the positive strand-
            # centric displacement needs to be reversed:
            return -posStrandDisp


def motifPosCmp(motifOcc1, motifOcc2):
    """Compares two motif occurences according to their position on a given
    sequence. The occurences must be on the same sequence."""
    assert(motifOcc1.getSeqName() == motifOcc2.getSeqName())

    motifOcc1Centre = motifOcc1.getCentre()
    motifOcc2Centre = motifOcc2.getCentre()

    if (motifOcc1Centre < motifOcc2Centre):
        return -1
    elif (motifOcc1Centre == motifOcc2Centre):
        return 0
    else:
        assert (motifOcc1Centre > motifOcc2Centre)
        return 1


def motifScoreCmp(motifOcc1, motifOcc2):
    """Compares two motif occurences according to their pval."""

    if (motifOcc1.getPval() < motifOcc2.getPval()):
        return -1
    elif (motifOcc1.getPval() == motifOcc2.getPval()):
        return 0
    else:
        assert (motifOcc1.getPval() > motifOcc2.getPval())
        return 1


def parse_MAST(mastFile, motifString = None):
    """Parse a MAST output file that is in text hitlist format, to produce a
    dictionary specifying all motif occurences."""

    # XXX Need to change this function interface, so that the function can
    # return sequence information (as well as match information) as well,
    # if requested. PROBLEM: That information is *not* contained in the
    # mast output; I am getting confused.

    motifOccurences = {}

    for line in mastFile.xreadlines():
        if (line[0] != "#"):
            elems = line.split()
            try:
                seqName = elems[0]
                strand = elems[1]
                start = int(elems[2])
                end = int(elems[3])
                pVal = float(elems[5]) # Changed on 5th Nov 2010: MAST changed!
                currHit = motifOccurence(seqName, strand, start, end, pVal, \
                                         motifName = motifString)
            except ValueError, e:
                print >> sys.stderr, "Invalid MAST output file line:", line
                raise e
            if (motifOccurences.has_key(seqName)):
                motifOccurences[seqName].append(currHit)
            else:
                motifOccurences[seqName] = [currHit]

    for seqName in motifOccurences.keys():
        motifOccurences[seqName].sort(motifPosCmp)

    return motifOccurences


def glam2scan2bed(motifFilename, sequenceFilename, nHits, outFile):
    """Runs GLAM2SCAN on the specified sequence, breaking the sequence
    up into smaller sections first, if it is too large. Output the motif hits
    in BED format, to the specified output file."""

    return # XXX FINISH IMPLEMENTING THIS FUNCTION LATER.
    # Calculate the required split sequence file size, based on memory usage
    # back-of-the-envelope calculation:
    # XXX

    # Break the input sequence into fragments that overlap by the motif length:
    # XXX

    # For each of the resulting sequence files...
    # XXX
    # Run glam2scan on the file:
    # XXX generate temporary glam2scan output filename
    # XXX run_GLAM2SCAN(sequenceFilename, motifFilename, outputFilename, 
    # nHits = nHits)
    #
    # Convert the output to bed format, printing to the specified output file:
    # XXX The BED (genomic) coordinates of a given site are calculated by adding the reported site position to the absolute genomic coordinates of the given sequence. Each resulting motif occurrence is appended as a line in the specified bed file.


def parse_GLAM2SCAN(glam2scanFile):
    """Parse a GLAM2SCAN output file, to produce a dictionary specifying all
    motif occurences."""

    motifOccurences = {}

    for line in glam2scanFile.xreadlines():
        elems = line.strip().split()
        if (len(elems) == 6):
            seqName = elems[0]
            # This line specifies a motif hit, => Store it:
            currHit = motifOccurence(seqName, elems[4], elems[1], elems[3], \
                                     None, None, float(elems[5]))
            if (motifOccurences.has_key(seqName)):
                motifOccurences[seqName].append(currHit)
            else:
                motifOccurences[seqName] = [currHit]

    for seqName in motifOccurences.keys():
        motifOccurences[seqName].sort(motifPosCmp)

    return motifOccurences


def writeFasta(seqs, outfile):
    """Writes the specified sequences in fasta format to the specified output
    file."""
    for seqName in seqs.keys():
        print >> outfile, ">" + seqName
        print >> outfile, seqs[seqName]


def parseFasta(seqFile):
    """Parses a fasta file, to produce a dictionary that has sequence names
    as keys and sequences as values."""

    seqs = {}
    currSeqname = None
    currSeq = ""
    for line in seqFile.xreadlines():
        if (line[0] == ">"):
            # The line designates a new sequence...

            # If there was a sequence previously, then add it to the dictionary
            # and start a new sequence out:
            if (currSeqname != None):
                seqs[currSeqname] = currSeq
                currSeq = ""
            elems = line.split()
            currSeqname = elems[0][1:]
        else:
            # The line is part of the current sequence => grow the current
            # sequence:
            currSeq = currSeq + line.strip()

    # Add the final sequence to the dictionary, if there was one:
    if (currSeqname != None):
        seqs[currSeqname] = currSeq

    return seqs


def parseFasta2(seqFile):
    """Parses a fasta file, to produce a dictionary that has sequence names
    as keys and arrays of sequences as values - for input datasets where
    multiple sequences can have the same name."""

    seqs = {}
    currSeqname = None
    currSeq = ""
    for line in seqFile.xreadlines():
        if (line[0] == ">"):
            # The line designates a new sequence...

            # If there was a sequence previously, then add it to the dictionary
            # and start a new sequence out:
            if (currSeqname != None):
                if (seqs.has_key(currSeqname)):
                    seqs[currSeqname].append(currSeq)
                else:
                    seqs[currSeqname] = [currSeq]
                currSeq = ""
            elems = line.split()
            currSeqname = elems[0][1:]
        else:
            # The line is part of the current sequence => grow the current
            # sequence:
            currSeq = currSeq + line.strip()

    # Add the final sequence to the dictionary, if there was one:
    if (currSeqname != None):
        if (seqs.has_key(currSeqname)):
            seqs[currSeqname].append(currSeq)
        else:
            seqs[currSeqname] = [currSeq]

    return seqs


# 30-09-2010: Created this function to wrap around pysam.faidx, to provide
# an easy and efficient way to retrieve genomic sequence data:
# 11-10-2010: Commenting out for the time being since I am not using it
# presently:
# def getSeqPysam(genomeFaFile, locBedItem):
#     """A wrapper function for calling pysam.faidx. genomeFaFile is the genome
#     in a single fasta file, and loc is a bed item specifying a chromosomal
#     location. Simply returns the returned value from pysam.faidx; i.e. a
#     list containing the fasta file output with one line per list element."""
#     # Convert the bed item into a chromosomal coordinate string:
#     # FIXME/PROBLEM: FaFile may have chromosomes numbered by number, rather than
#     # having "chr" at the start. How should I deal with this??
#     # For the time being, I will eliminate the "chr" at the start, thus assuming
#     # that the fasta file contains chromosome names as numbers/X/Y only:
#     chrom = locBedItem.getChrom()[3:]
#     chrStart = locBedItem.getChromStart()
#     chrEnd = locBedItem.getChromEnd()

#     queryString = chrom + ":" + str(chrStart) + "-" + str(chrEnd)

#     # Call pysam.faidx and process the resulting list, to generate a single
#     # sequence:
#     print >> sys.stderr, "Attempting to run pysam.faidx for the inputs", \
#         genomeFaFile, "and", queryString
#     return pysam.faidx(genomeFaFile, queryString)


def revcomp(input_seq):
    lett_converter = \
        {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', \
         'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'}
    
    tmp_seq = []
    for lett in input_seq:
        assert (lett in lett_converter)
        new_lett = lett_converter[lett]
        tmp_seq.append(new_lett)
        
    tmp_seq.reverse()
    output_seq = ""
    for lett in tmp_seq:
        output_seq = output_seq + lett

    return output_seq


def getSingleSequence(bedItem, database):
    """Retrieves a single sequence, with specified coordinates (as a BED_line
    object) from the specified database (e.g. mm8), and returns that
     sequence."""
    pid = os.getpgid(0)
    tmp_bed_filename = "bed_file_" + str(pid) + ".tmp"
    # Check if the file exists already:
    if (os.access(tmp_bed_filename, os.F_OK)):
        # The processed temporary filename is invalid, since a file of that
        # name already exists => Abort with an error:
        print >> sys.stderr, "ERROR: Proposed temporary file, " + tmp_bed_filename \
            + " already exists."
        sys.exit(1)
    tmp_bed_file = open(tmp_bed_filename, 'w')
    tmp_bed_file.write(bedItem.getChrom() + " " + str(bedItem.getChromStart()) \
                       + " " + str(bedItem.getChromEnd()) + " X\n")
    tmp_bed_file.flush()
    tmp_bed_file.close()


    # Create a temporary output filename:
    pid = os.getpgid(0)
    tmp_output_filename = "outputSeuqence_" + str(pid) + ".tmp"
    # Check if the file exists already:
    if (os.access(tmp_output_filename, os.F_OK)):
        # The processed temporary filename is invalid, since a file of that
        # name already exists => Abort with an error:
        print >> sys.stderr, "ERROR: Proposed temporary file, " + \
            tmp_output_filename + " already exists."
        sys.exit(1)

    # This function works by acting as a wrapper to getSequences:
    getSequences(tmp_bed_filename, database, tmp_output_filename)
    seq = parseFasta(open(tmp_output_filename))

    # Delete the temporary files that were used in this function:
    print >> sys.stderr, "Deleting temporary files,", \
       tmp_bed_filename, tmp_output_filename
    cmdresult = commands.getstatusoutput("rm " + tmp_bed_filename + " " + tmp_output_filename)
    print >> sys.stderr, cmdresult

    return seq


def getSequences(bedFileName, database, outputFileName):
  """Retrieves a fasta file of the sequences corresponding to genomic locations
  in the specified bed file. Retrieves the sequences from the specified
 database. A name for the output fasta file must be specified. This function
  causes that file to be filled with the retrieved sequences."""

  print >> sys.stderr, "Copying bed file to genome browser server. Command:"
  command = "scp " + bedFileName + " t.whitington@genome.imb.qfab.org:~/tmp"
  print >> sys.stderr, command
  cmdresult = commands.getstatusoutput(command)
  print >> sys.stderr, "Command Result: ", cmdresult

  print >> sys.stderr, "Retrieving sequences on genome browser server. Command:"
  command = "ssh t.whitington@genome.imb.qfab.org '/common/bin/x86_64/sequenceForBed -db=" + database + " -bedIn=tmp/" + bedFileName + " -fastaOut=tmp/" + outputFileName + "'"
  print >> sys.stderr, command
  cmdresult = commands.getstatusoutput(command)
  print >> sys.stderr, "Command Result: ", cmdresult

  print >> sys.stderr, "Copying resulting sequences from genome browser server. Command:"
  command = "rsync -v t.whitington@genome.imb.qfab.org:~/tmp/" + outputFileName + " ."
  print >> sys.stderr, command
  cmdresult = commands.getstatusoutput(command)
  print >> sys.stderr, "Command Result: ", cmdresult


# Replacing the old implementation of getBestFimoHits, to make it much faster,
# by using pandas:
# FIXME: NOTE: Fimo outfile is considered to be a filename now, not an open file
# object. May need to adapt calling code accordingly, although I actually think
# I have adjusted them all appropriately already.
def getBestFimoHits(fimoOutfilename, nHits=1, fimoType="fimo", chunkSize=5000000, tmpDir="/tmp"):
    """Parses the specified fimo outputfile and extracts the best hit for each
    sequence, returning a dictionary of sequenceName->motifName->[motifOccurence] objects."""

    # Process the file in blocks of around 10 million lines, to keep the memory
    # overhead low...

    seqName2motif2bestHits = {}

    # Skip the first line in the fimo outfile:
    fimoOutfile = open(fimoOutfilename)

    fimoOutfile.readline()

    # Code for breaking the fimo output file into chunks of around 10 million
    # lines and processing each chunk into the best hits dictionary. It's
    # ugly, but pretty sure this should work...

    currLine = fimoOutfile.readline()

    currChunkFilename = utility.makeTempFilename("currChunk_fimoData", tmpDir=tmpDir)

    totalLineCount = 0
    while currLine != "":
        # Extract and process current chunk...
        print >> sys.stderr, "Extracting current chunk of data. Chunksize:", chunkSize
        lineCounter = 0
        currChunkFile = open(currChunkFilename, 'w')
        while currLine != "" and lineCounter < int(chunkSize):
            print >> currChunkFile, currLine.strip()
            elems = currLine.strip().split()
            currLine = fimoOutfile.readline()
            lineCounter += 1
            totalLineCount += 1

        # Extract previous sequence and motif from current line:
        elems = currLine.strip().split()

        prevMotifName = None
        prevSeqName = None

        if currLine != "":
            prevMotifName = elems[0]
            prevSeqName = elems[1]

        currMotifName = prevMotifName
        currSeqName = prevSeqName

        # Keep reading lines until the current "block" of sequence match information
        # is completed:
        while currLine != "" and currMotifName == prevMotifName and currSeqName == prevSeqName:
            print >> currChunkFile, currLine.strip()
            currLine = fimoOutfile.readline()
            totalLineCount += 1
            elems = currLine.strip().split()

            currMotifName = elems[0]
            currSeqName = elems[1]
            
        # Process the extracted chunk of fimo data, updating the
        # seqName2motif2bestHits dictionary accordingly:
        print >> sys.stderr, "Total lines read so far:", totalLineCount
        print >> sys.stderr, "Processing current chunk..."
        currChunkFile.close()
        getBestHitsForChunk(currChunkFilename, seqName2motif2bestHits,
                            nHits=nHits, fimoType=fimoType)
        
    try:
        cmdlineProgs.deleteFiles([currChunkFilename])
    except Exception, e:
        pass # Don't worry if file didn't exist.

    return seqName2motif2bestHits


def getBestHitsForChunk(fimoChunkFilename, seqName2motif2bestHits, nHits=1, fimoType="fimo"):
    # Parse the fimo output file using pandas, and filter for
    # rows containing the best (lowest p-value) hit for each
    # motif/sequence pair...
    print >> sys.stderr, "Processing chunk:", fimoChunkFilename
    fimoDataFrame = pandas.read_csv(fimoChunkFilename, sep = "\t", header=None, names=["Motif", "Seq", "Start", "Stop", "Strand", "Log-odds", "p-value", "Site"])

    print >> sys.stderr, "Grouping lines..."
    grouped = fimoDataFrame.groupby(["Motif", "Seq"])

    # def idxOfMin(items):
    #     minVal = min(items)
    #     itemIdx = 0
    #     minItemIdxs = []
    #     for key in items.keys():
    #         if items[key] == minVal:
    #             minItemIdxs.append(str(key))
    #     itemIdx += 1
    #     return random.choice(minItemIdxs)

    print >> sys.stderr, "Getting best hit indexes..."
    #minPvalIdxs = grouped.agg({'p-value' : idxOfMin})
    minPvalIdxs = grouped.agg({'p-value' : np.argmin})

    print >> sys.stderr, "Getting best hits..."
    idxList = list(minPvalIdxs['p-value'])
    idxListInts = map(lambda val: int(val), idxList)
    bestHits = fimoDataFrame.iloc[idxListInts]

    # This step is slow. Surely there must be a way to improve on this!
    # Convert the above data frame of best hits into the required dictionary
    # structure...
    print >> sys.stderr, "Converting data frame into dictionary structure..."
    for idx in range(len(bestHits)):
        currHit = bestHits.iloc[idx]
        motif = currHit[0]
        seqName = currHit[1]
        bestHitAsLine = "%s\t%s\t%d\t%d\t%s\t%1.6g\t%1.6g\t%s" % tuple(currHit)
        currHitObj = parseFimoLine(bestHitAsLine, fimoType=fimoType)
        if not seqName2motif2bestHits.has_key(seqName):
            seqName2motif2bestHits[seqName] = {}

        seqName2motif2bestHits[seqName][motif] = [currHitObj]

    return seqName2motif2bestHits


# Introduced on April 5th 2013 for use in snp allele scanning with many motifs
# simultaneously.
def getBestFimoHitsOLD(fimoOutfile, nHits=1, fimoType="fimo"):
    """Parses the specified fimo outputfile and extracts the best hit for each
    sequence, returning a dictionary of sequenceName->motifName->[motifOccurence] objects."""

    # Modifying on June 19th 2013 to process it in blocks; one block per motif.
    # This was necessary to limit memory usage:
    seqName2motif2bestHits = {}
    currLine = fimoOutfile.readline()
    assert currLine[:5] == "Motif" # First line should be a header; discard
    currLine = fimoOutfile.readline()
    while currLine != "":
        (currLine, currSeq2pwmHit) = parseFimoBlock(currLine, fimoOutfile, fimoType)
        extractBestHits(currSeq2pwmHit, seqName2motif2bestHits, nHits=1)

    return seqName2motif2bestHits


def extractBestHits(seq2pwmHit, seqName2motif2bestHits, nHits=1):
    """Extracts the best hits from the specified block of fimo data, into
    the specified output dictionary."""

    # Use a hash of hashes to link each sequence to all hits for each
    # motif:
    seqName2motif2hits = {}

    # Process the hits into the above data structure:
    for seqName in seq2pwmHit.keys():
        seqName2motif2hits[seqName] = {}
        currSeqAllHits = seq2pwmHit[seqName]
        for currHit in currSeqAllHits:
            motifName = currHit.getMotifName()
            if seqName2motif2hits[seqName].has_key(motifName):
                seqName2motif2hits[seqName][motifName].append(currHit)
            else:
                seqName2motif2hits[seqName][motifName] = [currHit]

    # Transform that into the output data structure, which links
    # seqName to motif to the top nHits hits for this sequence/motif:
    for seqName in seqName2motif2hits.keys():
        if not seqName2motif2bestHits.has_key(seqName):
            seqName2motif2bestHits[seqName] = {}

        motifName = seqName2motif2hits[seqName].keys()[0]
        # Retrieve all hits for the current sequence/motif
        # combination:
        currHits = seqName2motif2hits[seqName][motifName]
        
        if (nHits == 1):
            # Get the single best hit; no sort, to make this
            # faster (otherwise I could just use the same
            # code as for (nHits > 1).
            
            # Get all equal best hits:
            currBestHits = [currHits[0]]
            hitIdx = 1
            while (hitIdx < len(currHits)):
                if (currHits[hitIdx].pVal < currBestHits[0].pVal):
                    currBestHits = [currHits[hitIdx]]
                elif (currHits[hitIdx].pVal == currBestHits[0].pVal):
                    # Hits are tied; keep all of them:
                    currBestHits.append(currHits[hitIdx])
                hitIdx = hitIdx + 1
            # Randomly choose one of the best hits from the one or more
            # equal-best hits identified:
            seqName2motif2bestHits[seqName][motifName] = [random.choice(currBestHits)]

        elif (nHits == -1):
            # Report all hits:
            currHits.sort(motifScoreCmp)
            seqName2motif2bestHits[seqName][motifName] = currHits
        else:
            # More general scenario; get the top nHits hits instead...
            
            # Sort the hits for the current sequence:
            currHits.sort(motifScoreCmp)
            seqName2motif2bestHits[seqName][motifName] = currHits[:nHits]


def getBestMASThits(mastOutfile, nHits=1):
    """Parses the specified mast outputfile and extracts the best hit for each
    sequence, returning a dictionary of sequenceName->motifOccurence objects."""

    allHits = parse_MAST(mastOutfile)
    bestHits = {}

    for seqName in allHits.keys():
        # Obtain the best hit from all hits for that sequence...
        currSeqAllHits = allHits[seqName]

        if (len(currSeqAllHits) == 0):
            bestHits[seqName] = None
        else:
            # Modifying on April 19th 2011: Generalising to multiple top hits...
            # Modifying on May 30th 2012: nHits of -1 means report all hits.

            if (nHits == 1):
                # Execute pre-existing code, to just get the single top hit...
                currBestHits = [currSeqAllHits[0]]
                hitIdx = 1
                while (hitIdx < len(currSeqAllHits)):
                    if (currSeqAllHits[hitIdx].pVal < currBestHits[0].pVal):
                        currBestHits = [currSeqAllHits[hitIdx]]
                    elif (currSeqAllHits[hitIdx].pVal == currBestHits[0].pVal):
                        # Hits are tied; keep all of them:
                        currBestHits.append(currSeqAllHits[hitIdx])
                    hitIdx = hitIdx + 1
                # Randomly choose one of the best hits from the one or more
                # equal-best hits identified:
                bestHits[seqName] = [random.choice(currBestHits)]
            elif (nHits == -1):
                # May 30th 2012: Report all hits.

                currSeqAllHits.sort(motifScoreCmp)
                bestHits[seqName] = currSeqAllHits
            else:
                # More general scenario; get the top nHits hits instead
                # (required as of April 19th 2011)...

                # Sort the hits for the current sequence:
                currSeqAllHits.sort(motifScoreCmp)
                bestHits[seqName] = currSeqAllHits[:nHits]
                
    return bestHits


def makeSeqAlnWiggleHtml(outFile, outDir, alignedSeqs, alignedWigProf,
                         filePrefix, ylim=None):
    """Inputs:
    - Output html filehandle.
    - Location of output directory.
    - A list of aligned sequences to generate a sequence logo from.
    - A wigProfile object representing the wiggle profile that is aligned with
    the sequences.

    Outputs:
    - Generates a file containing a sequence logo representing the aligned
    sequences.
    - A figure file showing the IC of the aligned sequences, aligned with the
    input wiggle profile.
    - Writes html to the specified file, drawing a little table containing the
    sequence logo and the IC-wiggle plot."""

    # Generate a sequence logo representing the sequence alignment...

    # Make a pwm from the alignment:
    tmpMotif = motifModule.pwm(alignedSeqs, dataType="seqAln")

    # Make the sequence logo:
    logoFilePrefix = outDir + "/" + filePrefix + "_logo"
    logoFileName = logoFilePrefix + ".png"
    tmpMotif.makeSeqLogo(logoFilePrefix, format="png")

    # Generate the IC-wiggle plot:
    icWiggleFilePrefix = outDir + "/" + filePrefix + "_icWiggle"
    icWiggleFilename = icWiggleFilePrefix + ".png"
    makeICwigglePlot(icWiggleFilePrefix, tmpMotif, alignedWigProf, ylim)

    # Generate html code making a table with elements linking to those two
    # files:
    logoAbsPath = os.path.normpath(\
        os.path.abspath(logoFileName))
    outDirAbsPath = os.path.normpath(os.path.abspath(\
            outDir))
    logoRelPath = utility.getRelPath(outDirAbsPath, logoAbsPath)

    icWiggleAbsPath = os.path.normpath(\
        os.path.abspath(icWiggleFilename))
    outDirAbsPath = os.path.normpath(os.path.abspath(\
            outDir))
    icWiggleRelPath = utility.getRelPath(outDirAbsPath, icWiggleAbsPath)

    table = '<table><tr><td><img src="./' + logoRelPath + \
        '" height="150"><br></td><td><img src="' + \
        icWiggleRelPath + \
        '" height="150"></tr></table>'
    outFile.write(table + '\n')


def makeICwigglePlot(filePrefix, motif, wigProf, ylim=None):
    """Generates an IC-wiggle plot, writing to a file with the specified prefix
    in its name"""

    # Generate barplot showing the specified wiggle profile:
    ax1 = pylab.subplot(212)
    ax1.cla()
    print >> sys.stderr, "Generating profile barplot..."
    profileSpan = range(wigProf.getProfLen())
    pylab.bar(profileSpan, wigProf.getProf(), color="blue")
#    pylab.xlim(-self.getFlankWidth(),
#                self.getFlankWidth()+motif.getWidth())
    if (ylim != None):
        # Use specified y limits:
        pylab.ylim(ylim[0], ylim[1])
    pylab.xlabel("Position (bp)", fontsize=20)
    pylab.ylabel("Avg wiggle score", fontsize=20)

    # Generate information content barplot to superimpose:
    ax2 = pylab.subplot(211)
    motifPosSpan = range(0, motif.getWidth())
    pylab.bar(motifPosSpan,
              motif.getIC_arr(), alpha=0.3,
              color="green")
    pylab.ylim(0,2)
#    pylab.xlim(-self.getFlankWidth(),
#          self.getFlankWidth()+motif.getWidth())
    pylab.ylabel("Motif column IC")
    
    profFilename = filePrefix + ".png"
    
    pylab.savefig(profFilename, format="png")
    ax1.cla()
    ax2.cla()
