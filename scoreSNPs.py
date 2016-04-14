#!/usr/bin/env python

import commands, math, pdb, sys
import genomics.formats.bed as bed
import getSeqFa, motifModule, seqUtils, utility, seqVariants
from optparse import OptionParser
import cmdlineProgs

class snvMotifEffect:
    def __init__(self, chrom, pos, snvName, motifName, refAlleleSeq,
                 bestAlleleSeq, worstAlleleSeq, bestHitStart,
                 bestHitEnd, bestHitStrand, bestHitSeq,
                 bestAlleleScore, worstAlleleScore, scoreDiff,
                 maxScore):
        self.chrom = chrom
        self.pos = pos
        self.snvName = snvName
        self.motifName = motifName
        self.refAlleleSeq = refAlleleSeq
        self.bestAlleleSeq = bestAlleleSeq
        self.worstAlleleSeq = worstAlleleSeq
        self.bestHitStart = bestHitStart
        self.bestHitEnd = bestHitEnd
        self.bestHitStrand = bestHitStrand
        self.bestHitSeq = bestHitSeq
        self.bestAlleleScore = bestAlleleScore
        self.worstAlleleScore = worstAlleleScore
        self.scoreDiff = scoreDiff
        self.maxScore = maxScore

        self.bestHitAvgCons = None
        self.matchedPeaks = []

    def setCons(self, avgCons):
        self.bestHitAvgCons = avgCons

    def findMatchedChipseq(self, nearbyChipseq, motifName2tfs, trace=False):
        """Given an array of nearby chip-seq peaks, finds all those
        that are for a TF corresponding to this motif, and stores them."""

        # Store the peaks information;
        # Convert the peaks list into a dictionary with tf name as key and
        # list of peaks as value:
        # FIXME: The chipseq peak items are just arrays; this is messy
        # and bug-prone.
        peaksDict = {}
        for peak in nearbyChipseq:
            tfName = peak[2]
            # TF called tfName is located near this SNV. Record this...
            if not peaksDict.has_key(tfName):
                peaksDict[tfName] = [peak]
            else:
                peaksDict[tfName].append(peak)

        # Get the corresponding TF names for this motif overlap, including
        # those joined by sequence similarity. This assumes that
        # self.motifName is a valid "motifPrefix" value:
        tfsForMotif = motifName2tfs[self.motifName]

        # See if any of the tfs (corresponding to this motif) are amongst the peaks
        # nearby this SNV:
        for tf in tfsForMotif:
            if peaksDict.has_key(tf):
                # The motif for this overlap matches the TF for a peak
                # close to the SNV:
                #peakInfo = reduce(lambda tok1, tok2: tok1 + "_" + tok2,
                #                  map(lambda tok: str(tok), peaksDict[tf]))

                peakInfo = peaksDict[tf]
                self.setMatched(peakInfo)

    def setMatched(self, peakInfo):
        # Hacking some changes on September 26th 2014; I think it should work
        # (for including more detailed info on the matched peaks):
        for peak in peakInfo:
            self.matchedPeaks.append(peak)

    def getBestAlleleScore(self):
        return self.bestAlleleScore

    def getMatched(self):
        return self.matchedPeaks

    def getMatchPeaksStr(self):
        # HACK: BROKEN: 24th October 2014:
        #ORIGINAL:
        #matchedPeakStrs = self.matchedPeaks
        matchedPeakStrs = map(lambda tok: str(tok), self.matchedPeaks[0])
        return reduce(lambda tok1, tok2: tok1 + "|" + tok2, matchedPeakStrs)

    def toString(self, outFormat="space", motifsCursor=None, snv=None, seqGetter=None):
        if outFormat == "space":
            effectStr = self.chrom + " " + str(self.pos) + " " + self.snvName + " " + self.motifName
        #effectStr = self.motifPrefix + " " + str(self.motifID) + " " + str(self.allele1minusLogP) + " " + str(self.allele2minusLogP) + " " + str(self.absDiff) + " " + str(self.chrom) + " " + str(self.pos) + " " + str(self.snvName) + " " + str(self.bestHitStart) + " " + str(self.bestHitEnd) + " " + str(self.avgCons) + " " + self.getMatchPeaksStr()
        elif (outFormat == "withMatchedPeaks"):
            effectStr = "(%s|%s)_%0.3f_%0.3f_%0.3f" % \
                (self.motifName, self.getMatchPeaksStr(), self.maxScore, self.scoreDiff, self.bestHitAvgCons)

        elif (outFormat == "withTFandGenomicSequence"):
            # Added on July 3rd 2014.
            # TFs_PWMname_bestAlleleReferenceSequencePlusStrand_bestAlleleScore

            # Need to be given access to the motifs database, the full
            # SNP object, and a seqGetter in order to do this. This is getting
            # quite ugly:
            assert (motifsCursor != None and snv != None and seqGetter != None)

            # Get list of corresponding TF HUGO names for the motif prefix:
            tfsForMotif = motifModule.getTFsForMotifPrefix(self.motifName,
                                                           motifsCursor)
            tfsString = "+".join(tfsForMotif)

            bestAlleleSeqMatchPlusStrand = self.getBestAlleleSeqMatchPlusStrand(snv, seqGetter)

            # HACK HERE TOO: Assuming that the score is a minus log 10 p-value.
            # Converting back to a p-value before printing:
            bestAlleleScore = self.getBestAlleleScore()
            assert(bestAlleleScore > 1)
            bestAllelePval = 10**(-bestAlleleScore)

            worstAlleleScore = self.worstAlleleScore
            #assert(worstAlleleScore > 1)
            worstAllelePval = 10**(-worstAlleleScore)
            
            effectStr = "%s#%s#%s#%s#%s#%1.3g#%s#%1.3g" % (tfsString, self.motifName,
                                                           bestAlleleSeqMatchPlusStrand,
                                                           self.bestHitStrand,
                                                           self.bestAlleleSeq,
                                                           bestAllelePval,
                                                           self.worstAlleleSeq,
                                                           worstAllelePval)
        else:
            assert outFormat == "underscore"
            effectStr = "(%s)_%0.3f_%0.3f_%0.3f" % \
                (self.motifName, self.maxScore, self.scoreDiff, self.bestHitAvgCons)
        return effectStr

    def getBestAlleleScore(self):
        return self.bestAlleleScore

    def getBestAlleleSeqMatchPlusStrand(self, snv, seqGetter):
        """Returns the genomic sequence corresponding to the best match of the
        given PWM to the best allele sequence."""

        # FIXME: Getting very hacky indeed. Needs a complete overhaul at some
        # stage (July 3rd 2014).

        # Figure out the index of the best allele within the snv:
        bestAllele = self.bestAlleleSeq
        alleleIdx = -1
        for alleleIdx in snv.getAlleleIdxs():
            if bestAllele == snv.getAlleleSeq(alleleIdx):
                bestAlleleIdx = alleleIdx

        # This might work but I'm not 100% sure (e.g. with indels):
        bestHitGenomicSequence = \
            seqGetter.getSequence(self.chrom + ":" + str(self.bestHitStart) + "-" + str(self.bestHitEnd))

        return bestHitGenomicSequence

    def toStringSmall(self):
        return self.motifPrefix + ";SCORES=" + str(self.allele1minusLogP) + "/" + str(self.allele2minusLogP) + "/" + str(self.absDiff) + ";CONS=" + str(self.bestHitAvgCons) + ";MATCHED=" + self.getMatchPeaksStr()

    def setConsFromReader(self, consReader):
        """Sets the conservation at this SNV motif overlap from the specifried
        conservation reader."""

        bestHitStart = self.bestHitStart
        bestHitEnd = self.bestHitEnd
                
        # Get the conservation data. Calculate the function (e.g.
        # average) over that data. If it passes the avg cons
        # threshold then record this hit:
        line = bed.BED_line(self.chrom + " " +
                            str(bestHitStart) +
                            " " + str(bestHitEnd))
        consData = consReader.getWig(line)
        # Hacking to deal with errors from inputs residing past the end of
        # chrosomes:
        if (len(consData) == 0):
            print >> sys.stderr, \
                "WARNING: Error with conservation data. Continuing anyway..."
            avgCons = 0
        else:
            avgCons = sum(consData)/float(len(consData))

        self.setCons(avgCons)


def parseSNVscores(scoresFile):
    """Parses the specified file (containing SNV motif score data) and returns
    a representation of this data"""

    snvMotifScores = {}

    for line in scoresFile.readlines():
        elems = line.strip().split("\t")
        posStr = elems[0]
        toks = posStr.split("_")
        chrom = toks[0]
        pos = int(toks[1])
        snvName = elems[1]
        motifName = elems[2]
        refAlleleSeq = elems[3]
        bestAlleleSeq = elems[4]
        worstAlleleSeq = elems[5]
        bestHitStart = int(elems[6])
        bestHitEnd = int(elems[7])
        bestHitStrand = elems[8]
        bestHitSeq = elems[9]
        bestAlleleScore = float(elems[10])
        worstAlleleScore = float(elems[11])
        scoreDiff = float(elems[12])
        maxScore = float(elems[13])

        currMotifEffect = \
            snvMotifEffect(chrom, pos, snvName, motifName, refAlleleSeq,
                           bestAlleleSeq, worstAlleleSeq, bestHitStart,
                           bestHitEnd, bestHitStrand, bestHitSeq,
                           bestAlleleScore, worstAlleleScore, scoreDiff,
                           maxScore)
        if snvMotifScores.has_key((chrom, pos)):
            snvMotifScores[(chrom, pos)].append(currMotifEffect)
        else:
            snvMotifScores[(chrom, pos)] = [currMotifEffect]

    return snvMotifScores


def scanSNVs(outFile, width2pwmFilenames, snvsFilename, seqGetter, options):
    # Scan all snvs in one go for each motif width...
    for width in width2pwmFilenames.keys():
        scanSNVsForWidth(outFile, snvsFilename, width, width2pwmFilenames[width], seqGetter,
                         scoreDiffThresh = options.scoreDiffThresh,
                         bestScoreThresh = options.bestScoreThresh,
                         showBestHitSeq = options.showBestHitSeq,
                         outputPvals = options.outputPvals,
                         tmpDir=options.tmpDir)


def groupMotifs(memeInfile, options, motifPseudo=0.01):
    width2pwmFiles = {}
    width2pwmFilenames = {}
    endOfFile = False

    while not endOfFile:
        try:
            currMotif = motifModule.pwm(memeInfile)
            currWidth = currMotif.getWidth()
            if (width2pwmFiles.has_key(currWidth)):
                outFile = width2pwmFiles[currWidth]
            else:
                tmpFilename = utility.makeTempFilename("Motifs_Width_" + str(currWidth), tmpDir=options.tmpDir)
                outFile = open(tmpFilename, 'w')
                width2pwmFiles[currWidth] = outFile
                width2pwmFilenames[currWidth] = tmpFilename
            currMotif.writeToMEME(outFile, pseudo=motifPseudo)
            print >> outFile, ""
        except ValueError, e:
            endOfFile = True
    for currWidth in width2pwmFiles.keys():
        width2pwmFiles[currWidth].close()

    return width2pwmFilenames


def tupCmp(idHitTup1, idHitTup2):
    """Specialised comparator function; tuples with two entries input. Compares
    based on the score of the second element in each."""
    if idHitTup1[1].getPval() < idHitTup2[1].getPval():
        return -1
    elif idHitTup1[1].getPval() == idHitTup2[1].getPval():
        return 0
    else:
        return 1


def scanSNVsForWidth(outFile, snvsFilename, width, memeFilename, seqGetter,
                     scoreDiffThresh = 0.5, bestScoreThresh = 3,
                     showBestHitSeq = False, outputPvals = False, tmpDir="/tmp"):
    """Scans all the snvs in the specified file with all the motifs in the
    specified meme file. All the motifs must be of the same specified width."""

    print >> sys.stderr, "Scanning SNVs for all motifs of width", width, "..."

    # Will store the SNV allele sequences too:
    snvAlleleRegSeqs = {}
    
    # Data structure to facilitate linking from snvName to SNV object:
    snvName2snv = {}

    # Start a snvAlleleSeq fasta file:
    snvAlleleSeqFilename = utility.makeTempFilename("snvAlleleSeqs_" + str(width),
                                                    check=False, tmpDir=tmpDir)
    snvAlleleSeqFile = open(snvAlleleSeqFilename, 'w')

    # Write out SNV allele sequences. Use iterative calls to
    # getAlleleWindowRegSeqs...
    #print >> sys.stderr, "Generating sequence file..."
    snvsFile = open(snvsFilename)
    for currLine in snvsFile.readlines():
        # Parse tokens into chrom, start, end etc. into a SNV object:
        elems = currLine.strip().split()
        chrom = elems[0]
        refPos = int(elems[1])
        name = elems[2]
        refSeq = elems[3]
        altSeqs = elems[4].split(",")
        alleles = [refSeq] + altSeqs
        currSNV = seqVariants.SNV(chrom, int(refPos), name, alleles)
        snvName2snv[currSNV.getLocStr()] = currSNV

        # Write all the SNV's allele sequences (given the current width)
        # to the snvAlleleSeq fasta file:
        currSNV.writeAlleleWindowSeqs(snvAlleleSeqFile, seqGetter, width)

        # FIXME: Copying/modifying blindly from processSNVs. This may not work.
        #for alleleIdx in currSNV.getAlleleIdxs():
        #    alleleID = currSNV.getAlleleID(alleleIdx)
        #    snvAlleleRegSeqs[alleleID] = alleleWindowRegSeqs[alleleIdx]

    # Flush all remaining lines to the snv allele sequences file:
    snvAlleleSeqFile.close()

    # Run fimo (with overlapping hits) with all the pwms on the SNV allele
    # sequence file. Get the best motif occurrence for all of the allele
    # sequences and pwms:
    print >> sys.stderr, "Running regular hacked fimo..."
    tmpFimo_outFilename = utility.makeTempFilename("fimoRun", check=False, tmpDir=tmpDir)
    print >> sys.stderr, snvAlleleSeqFilename
    seqUtils.runFimo(memeFilename, snvAlleleSeqFilename, tmpFimo_outFilename,
                     options = "--allow-overlap --output-pthresh 1",
                     verbose = True)

    print >> sys.stderr, "Processing fimo output..."
    alleleID2motif2bestHit = seqUtils.getBestFimoHits(tmpFimo_outFilename, tmpDir=tmpDir)
    print >> sys.stderr, "Done."

    # Cleanup: Delete the temporary SNV allele sequence file and fimo
    # output file:
    #cmdlineProgs.deleteFiles([snvAlleleSeqFilename, tmpFimo_outFilename])

    # Parse motifs here, in order to get motif names:
    memeInfile = open(memeFilename)
    motifsRemain = True
    motifs = []
    while motifsRemain:
        try:
            currMotif = motifModule.pwm(memeInfile)
            motifs.append(currMotif)
        except ValueError, e:
            motifsRemain = False

    # FIXME: Quite ugly code here. Some object-oriented approach might have
    # produced more modular code.

    # For each SNV/motif pair, get the best and worst allele hits,
    # compute/extract the relevant score statistics and write them to
    # the relevant output file...
    for snvName in snvName2snv.keys():
        #print >> sys.stderr, " TRACE: Outputting stats for SNV", snvName, "..."
        currSNV = snvName2snv[snvName]

        # Motif hits are not stored per SNV; get all hits for
        # this current SNV, organised by motifs.
        # For this SNV, for each motif, get all alleleIdx/hit tuples:
        motif2alleleHits = {}
        for alleleIdx in currSNV.getAlleleIdxs():
            alleleID = currSNV.getAlleleID(alleleIdx)
            for motif in motifs:
                motName = motif.getName()
                try:
                    alleleMotifBestHit = \
                        alleleID2motif2bestHit[alleleID][motName][0]
                except Exception, e:
                    pdb.set_trace()
                    dummy = 1
                if not motif2alleleHits.has_key(motName):
                    motif2alleleHits[motName] = \
                        [(alleleIdx, alleleMotifBestHit)]
                else:
                    motif2alleleHits[motName].append((alleleIdx,
                                                      alleleMotifBestHit))

        # Sort the alleleID/hit tuples according to hit score, for
        # each motif, so the best and worst can be retrieved:
        for motifName in motif2alleleHits.keys():
            motif2alleleHits[motifName].sort(tupCmp)

        # Output the relevant statistics for this SNV, for each motif...
        for motifName in motif2alleleHits.keys():
            #print >> sys.stderr, "  TRACE: Outputting stats for motif", motifName, "..."

            # Get the best allele's and worst allele's alleleIdx/motifHit
            # tuple:
            bestAlleleHitTup = motif2alleleHits[motifName][0]
            worstAlleleHitTup = motif2alleleHits[motifName][-1]

            # Get the allele sequences:
            bestAlleleIdx = bestAlleleHitTup[0]
            worstAlleleIdx = worstAlleleHitTup[0]
            bestAlleleSeq = currSNV.getAlleleSeq(bestAlleleHitTup[0])
            worstAlleleSeq = currSNV.getAlleleSeq(worstAlleleHitTup[0])

            # Get the motif hits for the best and worst alleles:
            bestAllele_hit = bestAlleleHitTup[1]
            worstAllele_hit = worstAlleleHitTup[1]

            # Hack: Calculate the position of the best allele's best binding
            # site here...
            
            # Calculate the start position of the scanned sequence:
            seqGenomicStart = currSNV.getRefPos() - width + 1
        
            # Calculate the genomic start and end position of the best binding site:
            seqStart = bestAllele_hit.getStart()
            seqEnd = bestAllele_hit.getEnd()
            bestHitGenomicStart = seqGenomicStart + seqStart - 1
            bestHitGenomicEnd = seqGenomicStart + seqEnd - 1

            # Potentially time-consuming due to disk access; get the best hit's
            # genomic sequence for this motif, if told to do so:
            bestHitCoordsStr = chrom + ":" + str(bestHitGenomicStart) + "-" + \
                str(bestHitGenomicEnd)
            bestAlleleHitSequence = None
            if (showBestHitSeq):
                bestAlleleFullSeq = \
                    currSNV.getAlleleWindowSeq(bestAlleleIdx, seqGetter, width)

                # October 2 2014: I think this was buggy; shouldn't have to take
                # reverse complement, as the variable bestAlleleFullSeq... Aha!
                # I should ?!? flip the PWM instead?? Now I'm confused...
                # XXX CONTINUE HERE.
                #if bestAllele_hit.getStrand() == "+":
                bestAlleleHitSequence = bestAlleleFullSeq[seqStart-1:seqEnd]
                #else:
                #    bestAlleleHitSequence = \
                #        seqUtils.revcomp(bestAlleleFullSeq[seqStart-1:seqEnd])

            bestAlleleScore = None
            worstAlleleScore = None
            scoreDiff = None
            maxScore = None

            if outputPvals:
                # Output scores as p-values and differences as p-value ratios,
                # instead of minus log values:
                bestAlleleScore = bestAllele_hit.getPval()
                worstAlleleScore = worstAllele_hit.getPval()
                scoreDiff = worstAlleleScore/bestAlleleScore
                bestScore = min(bestAlleleScore, worstAlleleScore)
            else:
                bestAlleleScore = -math.log(bestAllele_hit.getPval(), 10)
                worstAlleleScore = -math.log(worstAllele_hit.getPval(), 10)
                scoreDiff = bestAlleleScore - worstAlleleScore
                bestScore = max(bestAlleleScore, worstAlleleScore)
                
            if ((outputPvals and (scoreDiff >= float(scoreDiffThresh) and
                                  bestScore < float(bestScoreThresh))) or
                ((not outputPvals) and (scoreDiff >= float(scoreDiffThresh) and
                                        bestScore > float(bestScoreThresh)))):
                print >> outFile, "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%0.3g\t%0.3g\t%0.3g\t%0.3g" % (currSNV.getLocStr(), currSNV.getName(), motifName, currSNV.getRefAlleleSeq(), bestAlleleSeq, worstAlleleSeq, bestHitGenomicStart, bestHitGenomicEnd, bestAllele_hit.getStrand(), bestAlleleHitSequence, bestAlleleScore, worstAlleleScore, scoreDiff, bestScore)


def main():
    # Parse the command-line arguments...
    description = """usage: %prog [options] <pwmFile> <snvFile> <genomeFile>\n
Inputs:
- A MEME file containing one or more motifs, which will be scanned over the
specified SNVs
- A file specifying SNV information, in vcf format, with respect to the specified
genome file; The following fields are read:
-- Chrom
-- Reference genome coordinates of the first letter of the reference allele
-- Name of the SNV
-- Reference allele sequence
-- One or more alternative allele sequences, comma-separated
- A genome fa file from the same assembly as the SNV

Outputs:
- Tab-delimited output to stdout, containing the following information:
-- SNV chrom, start, end, name, strand, observed as in input
-- Motif name
-- Sequence of the reference allele
-- Sequence of preferred allele
-- Sequence of worst allele
-- Genomic start and end of the best hit
-- Genomic strand of best hit
-- Genomic sequence of the best hit if requested
-- highest -log10(p-value) of PWM scanned over best allele
-- highest -log10(p-value) of PWM scanned over worst allele
-- Difference between the above two scores
"""

    parser = OptionParser(usage = description)
#    parser.add_option("--mysqlUname", dest = "uname", default = "tom",
#                      help = "Mysql username for logging into the GOGE" + 
#                      "database. [default: %default]")
    parser.add_option("--motifPseudo", dest = "motifPseudo",
                      default = "0.01",
                      help = "Pseudo-count to add to motif column values " + \
                          "when parsing from meme file. Default=[%default]")
    parser.add_option("--outputPvals", action="store_true", dest="outputPvals",
                      default = False,
                      help = "Output motif scores as p-values and pval ratios.")
    parser.add_option("--tmpDir", dest="tmpDir", default = "/tmp",
                      help = "Location of the temporary directory, where tmp files will be stored." +
                      " [default: %default]")
    parser.add_option("--scoreDiffThresh", dest = "scoreDiffThresh",
                      default = "0.5",
                      help = "Difference in scores must exceed this in " + \
                          "order to be reported. Default=[%default]")
    parser.add_option("--bestScoreThresh", dest = "bestScoreThresh",
                      default = "3",
                      help = "Best score must exceed this in " + \
                          "order to be reported. Default=[%default]")
    parser.add_option("--showBestHitSeq", action="store_true",
                      dest="showBestHitSeq",
                      help = "Show the sequence of the best hit for each" + \
                          " SNV/motif pair.")
    parser.add_option("--debug", action="store_true", dest="debug",
                      help = "Debug the program using pdb.")
    (options, args) = parser.parse_args()

    # Parse the input parameters...

    if (options.debug):
        pdb.set_trace()

    # Make sure the 2 required input files exist:
    if (len(args) != 3):
        print >> sys.stderr, "WRONG # ARGS: ", len(args)
        parser.print_help()
        sys.exit(1)

    # Group all motifs into separate meme files by width; one file per unique
    # width:
    memeInfile = open(args[0])
    width2pwmFilenames=groupMotifs(memeInfile, options, motifPseudo=float(options.motifPseudo))

    # Get the snvs filename and sequence retriever:
    snvsFilename = args[1]
    seqGetter = getSeqFa.seqRetriever(args[2])

    outFile = sys.stdout

    # Print the header:
    #print >> outFile, "SNV_Loc\tSNV_Name\tMotif\tRefAllele\tBestAllele\tWorstAllele\tBestStart\tBestEnd\tBestStrand\tBestSeq\tBestAlleleScore\tWorstAlleleScore\tScoreDiff\tBestScore"

    # Scan the SNVs:
    scanSNVs(outFile, width2pwmFilenames, snvsFilename, seqGetter, options)

    # FIXME: This is hacky and potentially dangerous: Delete the pwm filenames
    # (which have now been used), as they should have been temporary files:
    for filename in width2pwmFilenames.values():
        cmdlineProgs.deleteFiles([filename])


if __name__ == "__main__":
    sys.exit(main())
