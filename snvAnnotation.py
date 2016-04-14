#!/usr/bin/env python

import commands, MySQLdb, os, pdb, re, sys, cPickle
import utility, motifModule, seqVariants, scoreSNPs, getSeqFa
import xml.etree.ElementTree as ET
import genomics.formats.bed as bed
import genomics.formats.wibIO as wibIO
import genomics.formats.wigIO as wigIO
import genomics.seqExpts.chipSeq as chipSeq
import chipseq2igvServerXML
import cmdlineProgs

class AnnotGroup:
    """Annotation of a group of SNVs with respect to a set of genomic
    features."""

    def __init__(self):
        # Links SNV object to SNV_Annotation object:
        self.snv2annot = {}

        # Links SNV location tuples (chrom, pos, name) to SNV object;
        # facilitates retrieval of annotation specified by SNV coords:
        self.loc2snv = {}

    def pickle(self, outfile):
        try:
            cPickle.dump(self, outfile)
        except Exception, e:
            print >> sys.stderr, e
            dummy = 1
            pdb.set_trace()

    def getAnnot(self, snv):
        return self.snv2annot[snv]

    def getAnnotFromLoc(self, snvLoc):
        try:
            return self.getAnnot(self.loc2snv[snvLoc])
        except Exception, e:
            print >> sys.stderr, "ERROR: NO LOC2SNV ENTRY:", snvLoc
            pdb.set_trace()
            x = 1

    def getAnnots(self):
        return self.snv2annot.values()

    def addSNVannot(self, snvAnnot):
        snv = snvAnnot.getSNV()
        self.snv2annot[snv] = snvAnnot
        self.loc2snv[snv.getLocTup()] = snv

    def getSNVs(self):
        return self.snv2annot.keys()

    def merge(self, annotGroup):
        """Merges the specified snv annotation object with this one."""

        for snv in annotGroup.getSNVs():
            assert not self.snv2annot.has_key(snv)
            currAnnot = annotGroup.getAnnot(snv)
            self.addSNVannot(currAnnot)

    def writeToCSV(self, outfile, motifsCursor, seqGetter, ldInfo = None):
        """Generates a text csv output table for this annotation.
Implements the following output format:

One row per analysed SNP

Columns:
* SNP name
* Comma-separated list of (tag, LD) tuples indicating which tag SNPs this SNP is in LD with
* List of DNase overlaps
* List of Footprint overlaps
* List of ChIP-seq nearby
* List of occupied ChIP-seq nearby
* List of occupied conserved ChIP-seq nearby
* List of bed feature overlaps
* DNase score and percentile tuples, comma-separated
* TF/best_PWM/bestHitReferenceSequencePlusStrand/max_score tuples, comma-separated
"""

        for snv in self.getSNVs():
            outputString = snv.toString()
            if ldInfo != None:
                tagSNVs = ldInfo.getTagPartners(snv)
                assert (len(tagSNVs) > 0)
                tagSNVsStr = reduce(lambda tok1, tok2: tok1 + "," + tok2,
                                    map(lambda tagSNV: tagSNV.toString(), tagSNVs))
                outputString = outputString + " " + tagSNVsStr

            outputString = outputString + " " + \
                self.getAnnot(snv).toCSVline(motifsCursor, seqGetter)
            print >> outfile, outputString

    def makeSNPlocBed(self, outFile):
        snpStringsOutput = {}
        for currSNV in self.getSNVs():
            bedString = currSNV.toBedString()
            if not snpStringsOutput.has_key(bedString):
                print >> outFile, bedString
            snpStringsOutput[bedString] = 1


class SNV_Annotator:
    """Annotates SNVs with respect to a set of genomic features."""

    def __init__(self, seqGetter, options):
        self.seqGetter = seqGetter

        self.options = options

        # Lower than lowest possible conservation score.
        # FIXME: Duplicated:
        self.lowCons = -10
    
        # Parse relevant options:
        self.checkMotifThresh = float(options.checkMotif)
        self.motifConsThresh = float(options.motifCons)
        self.matchPeakDist = int(options.matchPeak)
        self.peakDist = int(options.peak)
        self.chipseqDist = int(options.chipseqDist)
        self.snvConsThresh = float(options.snvCons)
    
        snpEffectDb = MySQLdb.connect(host=options.host,
                                      db="snpEffect",
                                      user=options.uname,
                                      read_default_file="~/.my.cnf")
        self.snpEffectCursor = snpEffectDb.cursor()

        motifsDb = MySQLdb.connect(host=options.host,
                                   db="motifs",
                                   user=options.uname,
                                   read_default_file="~/.my.cnf")

        self.motifsCursor = motifsDb.cursor()
        
        self.consReader = None
        if (self.motifConsThresh > self.lowCons or
            self.snvConsThresh > self.lowCons):
            self.consReader = wibIO.wibReader(options.consDir)

        
        # FIXME: Nasty, but keeping for the time being until I reimplement this whole
        # analysis more cleanly.
        # Used for getting integer chromosome code in mysql database tables.
        self.chrStr2chrInt = {"1":1, "2":2, "3":3, "4":4, "5":5, "6":6, "7":7, "8":8, "9":9, "10":10, "11":11, "12":12, "13":13, "14":14, "15":15, "16":16, "17":17, "18":18, "19":19, "20":20, "21":21, "22":22, "X":23, "Y":24, "chr1":1, "chr2":2, "chr3":3, "chr4":4, "chr5":5, "chr6":6, "chr7":7, "chr8":8, "chr9":9, "chr10":10, "chr11":11, "chr12":12, "chr13":13, "chr14":14, "chr15":15, "chr16":16, "chr17":17, "chr18":18, "chr19":19, "chr20":20, "chr21":21, "chr22":22, "chrX":23, "chrY":24}
    
        # Get all motifs to consider, in a data structure to facilitate effecient
        # SNP annotation...
        # Make a combined file of all motifs according to the specified
        # motif prefixes:
        # FIXME: Hard-coded
        motifBaseDir = options.motifBaseDir
        suffix = "_fw_ic0.meme"
        combinedMEMEfilename = utility.makeTempFilename("MEME_combined", tmpDir=self.options.tmpDir)
        combinedMEMEfile = open(combinedMEMEfilename, 'w')
        prefixes = map(lambda line: line.strip(),
                       open(self.options.motifPrefixes).readlines())
        motifModule.motifPrefixes2meme(combinedMEMEfile,
                                       prefixes,
                                       motifBaseDir,
                                       suffix)

        # Retrieve the set of TFs corresponding to each motif name, so that this
        # doesn't need to be done many times later:
        motifName2tfs = {}
        for motifName in prefixes:
            motifName2tfs[motifName] = \
                motifModule.getTFsForMotifPrefix(motifName,
                                                 self.motifsCursor)

        self.motifName2tfs = motifName2tfs

        combinedMEMEfile.close()
        self.width2pwmFilenames = scoreSNPs.groupMotifs(open(combinedMEMEfilename), self.options)
        cmdlineProgs.deleteFiles([combinedMEMEfilename])

        # June 29th: A bit hacky. Getting the specified set of motifs for
        # scanning SNPs that have DNase overlap:
        self.width2pwmFilenames_dnase = None
        if options.dnaseFocus != None and options.motifPrefixes_DNase != None:
            suffix = "_fw_ic0.meme"
            combinedMEMEfilename_dnase = utility.makeTempFilename("MEME_combined_DNase", tmpDir=self.options.tmpDir)
            combinedMEMEfile_dnase = open(combinedMEMEfilename_dnase, 'w')
            prefixes_dnase = map(lambda line: line.strip(),
                                 open(self.options.motifPrefixes_DNase).readlines())
            motifModule.motifPrefixes2meme(combinedMEMEfile_dnase,
                                           prefixes_dnase,
                                           motifBaseDir,
                                           suffix)
            combinedMEMEfile_dnase.close()
            self.width2pwmFilenames_dnase = scoreSNPs.groupMotifs(open(combinedMEMEfilename_dnase), self.options)
            cmdlineProgs.deleteFiles([combinedMEMEfilename_dnase])

        # Set up a data structure attribute, "bigWigs" recording the annotator's
        # bigWig files. Use a dictionary, linking from variableName to bigWig
        # filename:
        self.bigWigs = {}
        if options.bigWigs != None:
            for nameFilePairString in options.bigWigs.split(","):
                nameFilePair = nameFilePairString.split("=")
                name = nameFilePair[0]
                filename = nameFilePair[1]
                self.bigWigs[name] = filename

        # A bit hacky: Set up a second attribute, which will contain a map from
        # score to percentile for each bigWig variable:
        self.bigWigPercentiles = {}
        for variableName in self.bigWigs.keys():
            bigWigFilename = self.bigWigs[variableName]
            print >> sys.stderr, "Calculating percentiles for bigWig file", bigWigFilename
            self.bigWigPercentiles[variableName] = \
                wigIO.bigWig2percentiles(bigWigFilename)

    def findDNaseOverlap(self, currSNV):
        if self.options.regRegTissue != None:
            # Get overlap for the regulatory regions in the specified tissue:
            cmdStr = "select type, chrom, start, end, tissue from regRegion where MBRContains(chromPos, PointFromText('Point(" + str(self.chrStr2chrInt[currSNV.chrom]) + " " + str(currSNV.refPos) + ")')) and tissue like \"" + self.options.regRegTissue + "\";"
            self.snpEffectCursor.execute(cmdStr)
            rows = self.snpEffectCursor.fetchall()

            if self.options.trace == True:
                # Print out all rows of overlap:
                for row in rows:
                    print >> sys.stderr, "RegRegOverlap:", row

            regRegOverlaps = []
            if len(rows) > 0:
                regRegOverlaps.append(rows[0][4])
                regRegOverlaps = regRegOverlaps + map(lambda row: row[4], rows[1:])

            return regRegOverlaps
        else:
            return []

    def findBigWigAnnots(self, annotGroup):
        """Annotates the SNPs with the bigWig files for this annotator."""

        # Generate a temporary bed file of *all* SNPs in the specified
        # annotation group:
        snpsTmpBedFilename = utility.makeTempFilename("SNPs_tempBedFile", tmpDir=self.options.tmpDir)
        snpsTmpBedFile = open(snpsTmpBedFilename, 'w')
        annotGroup.makeSNPlocBed(snpsTmpBedFile)
        snpsTmpBedFile.close()

        for variableName in self.bigWigs.keys():
            # Get the bigWig filename:
            currBigWigFilename = self.bigWigs[variableName]

            # Generate a temporary output filename:
            wigScoresOutfilename = utility.makeTempFilename("wigScoresOutfile", tmpDir=self.options.tmpDir)

            dummyOutfilename = utility.makeTempFilename("dummyScoresOutfile", tmpDir=self.options.tmpDir)

            cmdlineProgs.runBigWigAverageOverBed(snpsTmpBedFilename,
                                                 currBigWigFilename,
                                                 wigScoresOutfilename,
                                                 dummyOutfilename,
                                                 verbose=True)

            # Parse the resulting output file, to generate a dictionary from
            # snpLocs to bigWig scores:
            currSNPloc2BigWigScore = wigIO.parseBigWigAvgOutput(open(wigScoresOutfilename))

            # Associate the snp locations' scores with the corresponding SNP
            # annotations:
            self.linkScoresToSNPannots(currSNPloc2BigWigScore, annotGroup, variableName)

            cmdlineProgs.deleteFiles([wigScoresOutfilename, dummyOutfilename])
        cmdlineProgs.deleteFiles([snpsTmpBedFilename])

    def linkScoresToSNPannots(self, snpLoc2BigWigScore, annotGroup, variableName):
        """Links up the scores to the specified snp annotations, for
        the specified variable name."""

        for currSNVloc in snpLoc2BigWigScore.keys():
            currSNVannot = annotGroup.getAnnotFromLoc(currSNVloc)
            currScore = snpLoc2BigWigScore[currSNVloc]
            try:
                currScorePercentile = self.bigWigPercentiles[variableName][currScore]
            except Exception, e:
                dummy = 1
                pdb.set_trace()

            currSNVannot.setScoreAndPercentile(variableName, currScore, currScorePercentile)

    # FIXME: This function was hacked together quickly to replace the previous
    # getAllChipseq() function, which used slow mysql range queries. It
    # is not elegant; e.g. it's confusing how I process both peaks and "bind poss"
    # at the same time. Some abstraction might help.
    def findAllChipseq(self, snv, peakDist):
        """Get all binding positions from ChIP-seq data. Use two tabix files,
        instead of using the mysql tables, since the tables are slow compared
        with tabix."""
    
        # ChIP-seq peaks...
    
        # Retrieve all chip-seq peaks nearby the snv:
        leftBorder = snv.refPos - peakDist
        rightBorder = snv.refPos + peakDist
    
        # Access the tabix file. Hacking to perform the fold-change/q-value
        # filter using awk. This was originally done using a "where" statement
        # in the corresponding mysql query:
        # FIXME: Should I pre-pend "chr" to the chromosome names or not??
        cmd = "tabix " + self.options.peaksTabixFile + " " + \
            str(snv.chrom) + ":" + str(snv.refPos - peakDist) + \
            "-" + str(snv.refPos + peakDist) + " | awk -F '#' '$(NF-1) < 5 && $NF > 2'"
        cmdResult = commands.getstatusoutput(cmd)
    
        outputLines = []
        if cmdResult[1] != '':
            outputLines = cmdResult[1].split("\n")
    
        linesAsToks = map(lambda outputLine: outputLine.replace("#", " ").split(), outputLines)
        chipseqInfo = []
        for toks in linesAsToks:
            currInfo = (toks[3], toks[4], toks[5], toks[0].replace("chr",""), int(toks[1]), float(toks[6]), float(toks[7]))
            chipseqInfo.append(currInfo)
    
        # Retrieve all chip-seq peaks nearby the snv:
        if (self.chipseqDist >= 0):
            # Won't implement this, since I won't use it yet anyway.
            # (WILL WHEN I SCAN PrCa REGIONS THOUGH!):
            pass
    
        return chipseqInfo
    
    def findFootprintOverlap(self, currSNV):
        if self.options.footprintTissue != None:
            # Get overlap for the footprints in this tissue:
            cmdStr = "select tissue from footprint where MBRContains(chromPos, PointFromText('Point(" + str(self.chrStr2chrInt[currSNV.chrom]) + " " + str(currSNV.refPos) + ")')) and tissue like \"" + self.options.regRegTissue + "\";"
            self.snpEffectCursor.execute(cmdStr)
            rows = self.snpEffectCursor.fetchall()

            if self.options.trace == True:
                # Print out all rows of overlap:
                for row in rows:
                    print >> sys.stderr, "Footprint:", row

            footprintOverlaps = []
            if len(rows) > 0:
                footprintOverlaps.append(rows[0][0])
                footprintOverlaps = footprintOverlaps + map(lambda row: row[0], rows[1:])

            return footprintOverlaps
        else:
            return []

    def filterSNVs(self, snvsFile, filteredSnvsFile, annotGroup,
                   matchCriteria = "chipseq", dnaseExpt = None):
        """Filter the SNVs in the specified snvs vcf file, retaining only those
        that have one or more features overlapping. Output to the specified output
        file. matchCriteria can be "chipseq" or "dnase"."""
        nOrig = 0 # Number of unfiltered SNVs
        nFiltered = 0 # Number of filtered SNVs
        for line in snvsFile.readlines():
            nOrig += 1
            elems = line.strip().split()
            chrom = elems[0]
            refPos = int(elems[1])
            name = elems[2]

            # Augst 5th. Filtering out weirdo snps that are named "esv[0-9]+":
            weirdo = False
            #re.match("esv[0-9]+", name) != None:
            if name[:3] == "esv":
                print >> sys.stderr, "Filtering out weirdo SNP:", name
                weirdo = True
    
            currAnnot = annotGroup.getAnnotFromLoc((chrom, refPos))
            
            if matchCriteria == "chipseq":
                # ChIP-seq filter:
                if ((len(currAnnot.getMatchedChipseq()) > 0) and (not weirdo)):
                    nFiltered += 1
                    print >> filteredSnvsFile, line.strip()
            else:
                # DNase filter:
                assert matchCriteria == "dnase"
                
                # Define the SNP as being overlapped by DNase if the specified
                # experiment is amongst the DNase elements overlapping this SNP,
                # or if the SNP coincides with DNase signal in the top 1% of
                # genomic signal levels for one of it's experiments:
                percentiles = map(lambda name: currAnnot.name2scoreAndPercentile[name][1],
                                  currAnnot.name2scoreAndPercentile.keys())
                snpOverlapped = \
                    (((len(currAnnot.getDNase()) > 0) and (dnaseExpt in currAnnot.getDNase())) \
                         or \
                         (max(percentiles) > 0.99))

                if (snpOverlapped and (not weirdo)):
                    nFiltered += 1
                    print >> filteredSnvsFile, line.strip()

        print >> sys.stderr, "Orig SNVs = %d. Filtered SNVs = %d." % \
            (nOrig, nFiltered)

    def scanSNVs(self, snvsFilename, annotGroup, width2pwmFilenames):
        """Scan the SNVs with the input motifs."""
    
        # Scan all filtered SNVs...
        print >> sys.stderr, "Scanning SNVs with motifs..."
        snvsScanFilename = utility.makeTempFilename("SNVmotifScan", tmpDir=self.options.tmpDir)
        snvsScanFile = open(snvsScanFilename, 'w')
        print >> sys.stderr, snvsScanFilename
        scoreSNPs.scanSNVs(snvsScanFile, width2pwmFilenames,
                           snvsFilename, self.seqGetter, self.options)
        snvsScanFile.close()
        snvsScanOutfile = open(snvsScanFilename)
        snvLoc2motifEffects = scoreSNPs.parseSNVscores(snvsScanOutfile)
        print >> sys.stderr, "Done."
        
        # HORRIBLE HACK: Added on July 28th 2014: Set the genomic sequence match
        # of the best-allele plus strand match, here:
        #for snvLoc in snvLoc2motifEffects:
        #    for currSnpMotifEffect in snvLoc2motifEffects[snvLoc]:
        #        currSNV = annotGroup.loc2snv[snvLoc]
        #        plusStrandMatch = \
        #            currSnpMotifEffect.getBestAlleleSeqMatchPlusStrand(currSNV, self.seqGetter)
        #        currSnpMotifEffect.bestHitSeq = plusStrandMatch

        # Cleanup; delete temporary width-grouped motif files:
        # FIXME: Is this working??
        cmdlineProgs.deleteFiles(width2pwmFilenames.values())

        # Delete filtered SNVs file:
        cmdlineProgs.deleteFiles([snvsScanFilename])

        return snvLoc2motifEffects

    def annotateSNV(self, snv):
        snvAnnot = SNV_Annotation(snv)
        
        # FIXME: Use VCF tools' intersect command instead? Might be much faster
        # and cleaner.
        
        # Get DNase and ChIP-seq overlap info...
        dnaseOverlap = self.findDNaseOverlap(snv)
        footprintOverlap = self.findFootprintOverlap(snv)
        
        # Only look for ChIP-seq overlap if told to:
        if self.options.peaksTabixFile != None:
            chipseqOverlap = self.findAllChipseq(snv, self.peakDist)
            matchChipseqOverlap = self.findAllChipseq(snv,
                                                      self.matchPeakDist)
        else:
            chipseqOverlap = []
            matchChipseqOverlap = []

        snvAnnot.setDNase(dnaseOverlap)
        snvAnnot.setFootprint(footprintOverlap)
        snvAnnot.setChipseq(chipseqOverlap)
        snvAnnot.setMatchedChipseq(matchChipseqOverlap)

        return snvAnnot

    def findBEDoverlaps(self, snvsFilename, annotGroup, bedFilename):
        """Annotates all specified SNVs (currently specified twice, in SNVs
        filename and annotGroup) with respect to the specified bed file."""

        # Run bedtools intersect on the SNVs file and the bed file:
        outFilename = utility.makeTempFilename("BedtoolsIntersectOutput", tmpDir=options.tmpDir)
        cmdlineProgs.runBedtoolsIntersect(bedFilename, snvsFilename, outFilename,
                                          verbose=True)

        # Parse the intersection output file and annotate the snvs accordingly:
        for line in open(outFilename).readlines():
          elems = line.strip().split()
          
          # Get chrom, start, end, and feature name:
          chrom=elems[0]
          start=elems[1]
          end=elems[2]
          featureName=elems[3]

          # Get the corresponding snv annot object from the annotGroup input
          # object...
          # HACK: This step is prone to error when the SNP directly coincides
          # with the exact end of the exon. In those cases, use "start"
          # as the snp position, rather than "start + 1". Need to engineer
          # this better:
          overlappedAnnot = None
          if annotGroup.loc2snv.has_key((chrom, int(start)+1)):
              overlappedAnnot = annotGroup.getAnnotFromLoc((chrom, int(start)+1))
          elif annotGroup.loc2snv.has_key((chrom, int(start))):
              overlappedAnnot = annotGroup.getAnnotFromLoc((chrom, int(start)))

          # Only register this annotation if the SNP was successfully retrieved
          # above; discard weirdo cases (e.g. where the SNP is a huge
          # "esv" variant):
          if overlappedAnnot != None:
              # Register this bed overlap for this annotation:
              overlappedAnnot.addBEDannot(bedFilename + "#" + featureName)

        cmdlineProgs.deleteFiles([outFilename])


    def annotateSNVs(self, snvsFilename):
        """Annotate the specified SNVs with respect to all features specified
        in the options. Returns an AnnotGroup object representing the
        annotations of all those SNVs."""

        # The resulting annotation:
        annotGroup = AnnotGroup()
    
        # For each SNV...
        print >> sys.stderr, "Annotating SNVs with DNase and ChIP-seq info..."
        snvsFile = open(snvsFilename)
        snvNum = 0
        for line in snvsFile.readlines():
            if snvNum % 100 == 0:
                print >> sys.stderr, "Processed %d SNPs." % snvNum
            snvNum += 1

            elems = line.strip().split()
            chrom = elems[0]
            refPos = elems[1]
            name = elems[2]
            alleleSeqs = [elems[3]] + elems[4].split(",")
            currSNV = seqVariants.SNV(chrom, int(refPos), name, alleleSeqs)

            currSNVannot = self.annotateSNV(currSNV)
            annotGroup.addSNVannot(currSNVannot)

        print >> sys.stderr, "Done."

        # Annotate all SNVs with the specified additional bed files (if any),
        # all in one go...

        # Parse the list of zero or more feature bed files:
        print >> sys.stderr, "Annotating all SNVs with bed features..."
        if self.options.bedFeatures != "NA":
            bedFilenames = self.options.bedFeatures.split(",")
            for bedFilename in bedFilenames:
                self.findBEDoverlaps(snvsFilename, annotGroup, bedFilename)
            print >> sys.stderr, "Done."


        # Do bigWig annotation...
        if self.options.bigWigs != None:
            # Annotate all SNPs with the bigWigs
            self.findBigWigAnnots(annotGroup)


        # Updated on June 29th to do ChIP-seq-overlapped scanning and
        # DNase-overlapped scanning separately.

    
        # ChIP-seq overlapped SNP scanning...
    
        # Filter snvs for chipseq overlap prior to running motif scan...
        print >> sys.stderr, "Filtering the SNVs prior to motif scanning..."
        filteredSnvsFilename = utility.makeTempFilename("FilteredSNVs", tmpDir=self.options.tmpDir)
    
        filteredSnvsFile = open(filteredSnvsFilename, 'w')
        snvsFile = open(snvsFilename)
    
        # Do the actual filter:
        self.filterSNVs(snvsFile, filteredSnvsFile, annotGroup,
                        matchCriteria = "chipseq")
        filteredSnvsFile.close()

        snvLoc2motifEffects = self.scanSNVs(filteredSnvsFilename, annotGroup,
                                            self.width2pwmFilenames)

        cmdlineProgs.deleteFiles([filteredSnvsFilename])
    
        # DNase (or other regulatory feature) overlapped SNP scanning...

        # Filter snvs for chipseq overlap prior to running motif scan...
        print >> sys.stderr, "Filtering the SNVs prior to motif scanning..."
        filteredSnvsFilename = utility.makeTempFilename("FilteredSNVs", tmpDir=self.options.tmpDir)
    
        filteredSnvsFile = open(filteredSnvsFilename, 'w')
        snvsFile = open(snvsFilename)
    
        # Do the actual filter:
        if self.options.dnaseFocus != None:
            self.filterSNVs(snvsFile, filteredSnvsFile, annotGroup,
                            matchCriteria = "dnase", dnaseExpt = self.options.dnaseFocus)
            filteredSnvsFile.close()

            snvLoc2motifEffects_dnase = self.scanSNVs(filteredSnvsFilename, annotGroup,
                                                      self.width2pwmFilenames_dnase)

            self.linkMotifsToSNPannots_dnase(snvLoc2motifEffects_dnase, annotGroup)

        cmdlineProgs.deleteFiles([filteredSnvsFilename])

        # Annotate each of the ChIP-seq-overlapped SNP motif matches...
    
        # Says whether to throw away motifEffect items that don't match
        # a corresponding TF's ChIP-seq. FIXME: Presently hard-coded:
        filterMotifEffects = True

        # Get conservation information and TF ChIP-seq match info for each of
        # the snv motif effects:
        # NOTE/FIXME: This might involve many iterations (due to many motif
        # hits), and could be a bottle-neck in the program.
        print >> sys.stderr, "Getting conservation and TF ChIP-seq match info..."
        for currSNVloc in snvLoc2motifEffects.keys():
            currSNVannot = annotGroup.getAnnotFromLoc(currSNVloc)

            motifEffectsArr = snvLoc2motifEffects[currSNVloc]
            newMotifEffectsArr = [] # Facilitates filtering
            for motifEffect in motifEffectsArr:
                # Told to score motif conservation?:
                if (self.motifConsThresh > self.lowCons):
                    motifEffect.setConsFromReader(self.consReader)
                # Told to match motif effects to TF ChIP-seq peaks?:
                if (self.matchPeakDist >= 0):
                    # Get all chip-seq peaks near the current snv:
                    nearbyChipseq = currSNVannot.getMatchedChipseq()
                    motifEffect.findMatchedChipseq(nearbyChipseq, self.motifName2tfs,
                                                   trace=self.options.trace)
                    # If told to do so, discard any motif effect
                    # that does not have a matched chipseq:
                    if (filterMotifEffects):
                        if len(motifEffect.getMatched()) > 0:
                            newMotifEffectsArr.append(motifEffect)

            if filterMotifEffects:
                if self.options.trace:
                    print >> sys.stderr, \
                        "Updating motif effects array; orig length = %d, filtered length = %d." % \
                        (len(snvLoc2motifEffects[currSNVloc]),
                         len(newMotifEffectsArr))

                # Have filtered out non-matched motifEffect items. Update
                # records accordingly:
                currSNVannot.setMotifEffects(newMotifEffectsArr)
            else:
                currSNVannot.setMotifEffects([])
        print >> sys.stderr, "Done."

        return annotGroup

    # HACKY (June 29th 2014):
    def linkMotifsToSNPannots_dnase(self, snvLoc2motifEffects, annotGroup):
        """For each of the motif effects, stored by snp loc, assigns them to the
        corresponding SNP annotation."""

        for currSNVloc in snvLoc2motifEffects.keys():
            currSNVannot = annotGroup.getAnnotFromLoc(currSNVloc)

            motifEffectsArr = snvLoc2motifEffects[currSNVloc]
            currSNVannot.setMotifEffects_dnase(motifEffectsArr)


# What a horrible mess: Introduced on July 28th 2014 to cope with JSON
# conversion for use in JavaScript (see SNV_Annotation.toAnnotForJSON()). However, it
# doesn't really make sense to implement something nice and neat until I am
# better at using javascript.
class SNV_AnnotForJSON:
    def __init__(self, snvAnnot):
        self.snv = snvAnnot.snv
        self.annotFeatures = AnnotFeatures(snvAnnot)


class AnnotFeatures:
    def __init__(self, snvAnnot):
        self.dnase = snvAnnot.dnase
        self.footprint = snvAnnot.footprint
        self.chipseq = snvAnnot.chipseq
        #self.matchChipseq = snvAnnot.matchChipseq
        self.motifEffects = snvAnnot.motifEffects
        self.motifEffects_dnase = snvAnnot.motifEffects_dnase
        self.bedOverlaps = snvAnnot.bedOverlaps
        self.name2scoreAndPercentile = snvAnnot.name2scoreAndPercentile


class SNV_Annotation:
    def __init__(self, snv):
        self.snv = snv

        # DNase regions overlapping it:
        self.dnase = []
    
        # DNase footprints overlapping it:
        self.footprint = []
    
        # Array of chipseqs nearby it:
        self.chipseq = []
    
        # ChIP-seq peak mappings to be used in motif + chips-seq
        # pairing:
        self.matchChipseq = []
    
        # Motif effects for the given SNV:
        self.motifEffects = []

        # Motif effects for the SNP, will only be set if it overlaps
        # DNase of a specified cell-type:
        self.motifEffects_dnase = []

        self.bedOverlaps = []

        # Flexible dictionary, linking variableName to variable score.
        # Introducing on June 30th 2014, to hold the scores from bigWig
        # files. But could also be used for other name->score mappings,
        # I suppose (?):
        self.name2scoreAndPercentile = {}

    def toAnnotForJSON(self):
        """Returns an alternative representation of this object. This is
        necessary because I want to have a single data structure containing
        all features, with feature name as key and feature value as value,
        rather than having each feature as a separate field of the SNP
        annotation (which would be produced by calling json.dump).

        This is a mess, and should be refactored at some point."""

        return SNV_AnnotForJSON(self)

    def setScoreAndPercentile(self, name, score, percentile):
        """Scores and percentiles for variable names. Introduced to deal
        with annotations from bigWig data."""

        self.name2scoreAndPercentile[name] = (score, percentile)

    def getBestScorePercentile(self):
        highestPercentile = -1
        for name in self.name2scoreAndPercentile.keys():
            scoreAndPercentile = self.name2scoreAndPercentile[name]
            if scoreAndPercentile[1] > highestPercentile:
                highestPercentile = scoreAndPercentile[1]

        return highestPercentile

    def getScoreAndPercentile(self, name):
        return self.name2scoreAndPercentile[name]

    def setDNase(self, dnaseOverlap):
        self.dnase = dnaseOverlap

    def getDNase(self):
        return self.dnase

    def getFootprint(self):
        return self.footprint

    def getChipseq(self):
        return self.chipseq

    def setFootprint(self, footprintOverlap):
        self.footprint = footprintOverlap

    def setChipseq(self, chipseqOverlap):
        self.chipseq = chipseqOverlap

    def setMatchedChipseq(self, matchChipseqOverlap):
        self.matchChipseq = matchChipseqOverlap

    def setMotifEffects(self, motifEffects):
        self.motifEffects = motifEffects

    def setMotifEffects_dnase(self, motifEffects):
        self.motifEffects_dnase = motifEffects

    def addBEDannot(self, bedFeature):
        """Registers a bed feature annotation for this SNP. Currently,
        the input argument is just a string.."""

        self.bedOverlaps.append(bedFeature)

    def getBEDannots(self):
        return self.bedOverlaps

    def getMatchedChipseq(self):
        return self.matchChipseq

    def getSNV(self):
        return self.snv

    def getBestMotifEffect(self):
        """Retrieves the motif effect with the highest best-allele score."""

        if len(self.motifEffects) > 0:
            bestMotifEffect = None
            bestMotifEffectScore = -1
            motEffectIdx = 0
            while motEffectIdx < len(self.motifEffects):
                currMotEffect = self.motifEffects[motEffectIdx]
                currMotEffectScore = currMotEffect.getBestAlleleScore()
                if currMotEffect.getBestAlleleScore() > bestMotifEffectScore:
                    bestMotifEffect = currMotEffect
                    bestMotifEffectScore = currMotEffectScore
                motEffectIdx += 1
            return bestMotifEffect
        else:
            return None

    def toCSVline(self, motifsCursor, seqGetter):
        """Generates a space-separated line of text, summarising all of the features
        for this annotation."""

        if len(self.dnase) > 0:
            dnaseOverlapsStr = reduce(lambda str1, str2: str1 + "," + str2,
                                      self.dnase)
        else:
            dnaseOverlapsStr = ""

        if len(self.footprint) > 0:
            footprintOverlapsStr = reduce(lambda str1, str2: str1 + "," + str2,
                                          self.footprint)
        else:
            footprintOverlapsStr = ""

        chipseqOverlapsStr = str(len(self.chipseq))

        motifEffectsOverlapsStrArr = \
            map(lambda sme: sme.toString(outFormat="withMatchedPeaks", motifsCursor = motifsCursor, snv = self.snv, seqGetter = seqGetter), self.motifEffects)

        if len(motifEffectsOverlapsStrArr) > 0:
            motifEffectsOverlapsStr = reduce(lambda str1, str2: str1 + "," + str2,
                                             motifEffectsOverlapsStrArr)
        else:
            motifEffectsOverlapsStr = ""

        # Modified on November 7th: Also calculate the most significant score
        # amongst *all* motif effects:
        motifEffectsBestAlleleScores = \
            map(lambda sme: sme.maxScore, self.motifEffects)
        maxScore = -1
        if len(motifEffectsBestAlleleScores) > 0:
            maxScore = max(motifEffectsBestAlleleScores)

        if (len(self.bedOverlaps) > 0):
            bedOverlapsStr = reduce(lambda str1, str2: str1 + "," + str2,
                                    self.bedOverlaps)
        else:
            bedOverlapsStr = ""

        # July 3rd 2014: Generate output for DNase score; name/score/percentile
        # tuples, comma-separated, but only if overlapped by relevant experiment
        # or if in top two percent:
        dnaseScoreOutputTups = []
        for dnaseName in self.name2scoreAndPercentile.keys():
            scoreAndPercentile = self.name2scoreAndPercentile[dnaseName]
            score = scoreAndPercentile[0]
            percentile = scoreAndPercentile[1]
            currTupString = "%s/%d/%1.4f" % (dnaseName, score, percentile)
            if percentile > 0.98 or len(self.dnase) > 0:
                dnaseScoreOutputTups.append(currTupString)
        dnaseScoreOutputStr = ",".join(dnaseScoreOutputTups)

        # July 3rd 2014: Generate output for affected PWM hits in DNase regions;
        # TF/best_PWM/bestHitReferenceSequencePlusStrand/max_score tuples,
        # comma-separated:
        motifEffectsDNaseOutputTups = map(lambda sme: sme.toString(outFormat="withTFandGenomicSequence",
                                                                   motifsCursor=motifsCursor,
                                                                   snv=self.snv,
                                                                   seqGetter=seqGetter),
                                          self.motifEffects_dnase)

        motifEffectsDNaseOutputStr = ",".join(motifEffectsDNaseOutputTups)

        motifEffectsDNaseBestAlleleScores = \
            map(lambda sme: sme.maxScore, self.motifEffects_dnase)
        maxScore_DNase = -1
        if len(motifEffectsDNaseBestAlleleScores) > 0:
            maxScore_DNase = max(motifEffectsDNaseBestAlleleScores)

        return dnaseOverlapsStr + " " + footprintOverlapsStr + " " + \
            chipseqOverlapsStr + " " + motifEffectsOverlapsStr + " " + \
            bedOverlapsStr + " " + str(maxScore) + " " + dnaseScoreOutputStr + \
            " " + motifEffectsDNaseOutputStr + " " + str(maxScore_DNase)

    def generateIGVbatchCmds(self, batchCmdsFile, igvSessionFilename, baseDir, outDir, height=1150):
        """Generates IGV batch commands text file contents for this
        annotation."""

        print >> batchCmdsFile, "new"
        print >> batchCmdsFile, "maxPanelHeight " + str(height)
        print >> batchCmdsFile, "load " + baseDir + "/" + igvSessionFilename
        print >> batchCmdsFile, "genome hg19"
        print >> batchCmdsFile, "snapshotDirectory " + baseDir + "/" + "/" + outDir
        annotSNV = self.getSNV()
        print >> batchCmdsFile, "region " + annotSNV.getChrom() + ":" + str(annotSNV.getRefPos()) + "-" + str(annotSNV.getRefPos())
        print >> batchCmdsFile, "goto " + annotSNV.getChrom() + ":" + str(annotSNV.getRefPos() - 500) + "-" + str(annotSNV.getRefPos() + 501)
        crmViewFilename = annotSNV.getChrom() + "_1kb-view_" + str(annotSNV.getRefPos()) + ".png"
        print >> batchCmdsFile, "snapshot " + crmViewFilename
        dnaViewFilename = annotSNV.getChrom() + "_50bp-view_" + str(annotSNV.getRefPos()) + ".png"
        print >> batchCmdsFile, "snapshot " + dnaViewFilename
        print >> batchCmdsFile, "goto " + annotSNV.getChrom() + ":" + str(annotSNV.getRefPos() - 25) + "-" + str(annotSNV.getRefPos() + 26)
        print >> batchCmdsFile, "snapshot " + snapshotName
        return (crmViewFilename, dnaViewFilename)

    # An evolving mess (July 25th 2014):
    def addFeaturesForIGV(self, featuresDict, featuresList, baseDataDir2urlBase):
        """Generates an IGV session file for visualising this snv annotation."""

        # Get all chipseq peak-calling analyses corresponding to the experiments
        # that have peaks nearby the SNV...

        # Get big nasty ugly tuples that include the identity of the chip-seq
        # experiment, as a string:
        nearbyPeaksTuples = self.getChipseq()
        nearbyPeaksExptStrings = map(lambda peakTup: peakTup[0].replace("|", "/"), nearbyPeaksTuples)

        # Load up the peak-calling analyses for those experiments:
        # FIXME: BROKEN: Currently, using a hack to determine chip-seq base
        # directory from expt string:
        
        def exptString2peakCallingLoc(exptStr):
            exptStrStart = exptStr.split("/")[0]
            # Hack/attempt to fix on September 1st:
            chipSeqBaseDir=baseDataDir2urlBase[1]#options.chipseqBaseDir
            if exptStrStart in ("HAIB", "SYDH", "UW"):
                chipSeqBaseDir="/home/tom/Work/data/processed/hg19/chipseq/ENCODE/"
            return chipSeqBaseDir + "/" + exptStr + "/Analyser.pk"

        nearbyPeaksPeakCallingAnalysisLocs = set(map(exptString2peakCallingLoc, nearbyPeaksExptStrings))
        peakCallingAnalysers = map(lambda analysisLoc: cPickle.load(open(analysisLoc)), nearbyPeaksPeakCallingAnalysisLocs)

        # For each of those peak-callings, generate a resource node and an xml node for loading that
        # peak-calling's relevant tracks:

        peakCallingAnalyses = map(lambda analyser: analyser.getAnalyses().values()[0], peakCallingAnalysers)
        peakCallingAnalyses.sort(chipSeq.cmpAnalyses)

        for peakCallingAnalysis in peakCallingAnalyses:
            # Specify the peaks as a resource:
            exptName = peakCallingAnalysis.getExpt().getName()
            #peakCallingAnalysis.writePeakCallsIGV_XML(resourcesNode, baseDataDir2urlBase, exptName)

            # Attribute value pairs for generating profile resource XML elements:
            # FIXME: Prone to being broken:
            profParamsTups = [(wigIO.ProfileParams.combined, 100)]
            att2val = chipseq2igvServerXML.getAtt2Val(peakCallingAnalysis, profParamsTups)

            # Specify the "combined" profile as a resource:
            profParams = wigIO.ProfileParams(100, wigIO.ProfileParams.combined)
            combinedProf = peakCallingAnalysis.getProfile(profParams)

            tdfLocAbs = os.path.abspath(combinedProf.getTDFLoc())
            # Horrible hack: April 13th 2015: Need to change URL for two distinct tdf location listings!:
            # Another compensatory horrible hack, April 30th 2015:
            urlLoc = tdfLocAbs
            if not tdfLocAbs.split('/')[9] == "SYDH":
                urlLoc1 = tdfLocAbs.replace(baseDataDir2urlBase[0], baseDataDir2urlBase[1])
                urlLoc = urlLoc1.replace("/home/bcbio/tmpMount", "/mnt/skjoldNoBackup")

            featuresDict[urlLoc] = {"name" : exptName}
            featuresList.append(urlLoc)

            for att in att2val:
                featuresDict[urlLoc][att] = att2val[att]
            
            #combinedProf.writeIGV_XML(resourcesNode, baseDataDir2urlBase, exptName, att2val)

        return len(peakCallingAnalyses)


class SNV_LDinfo:
    """Contains LD info for a particular set of tag SNVs and LD partner SNVs."""

    def __init__(self, ldInfoFile):
        # These map from a (chrom, pos) tuple to an array of (chrom, pos)
        # tuples:
        self.tagTup2ldPartnerTups = {}
        self.ldPartnerTup2tagTups = {}

        for line in ldInfoFile.readlines():
            elems = line.strip().split()
            tagSNVtup = (elems[0], int(elems[1]), elems[2])
            rSquareStr = elems[3]
            rSquare = -1
            if (rSquareStr != None):
                rSquare = float(rSquareStr)
            ldPartnerSNVtup = (elems[4], int(elems[5]), elems[6])

            if self.ldPartnerTup2tagTups.has_key(ldPartnerSNVtup):
               self.ldPartnerTup2tagTups[ldPartnerSNVtup].append((tagSNVtup, rSquare))
            else:
                self.ldPartnerTup2tagTups[ldPartnerSNVtup] = [(tagSNVtup, rSquare)]

            if self.tagTup2ldPartnerTups.has_key(tagSNVtup):
                self.tagTup2ldPartnerTups[tagSNVtup].append((ldPartnerSNVtup, rSquare))
            else:
                self.tagTup2ldPartnerTups[tagSNVtup] = [(ldPartnerSNVtup, rSquare)]

    def getLDpartnerTups(self, tagSNVtup):
        return self.tagTup2ldPartnerTups[tagSNVtup]

    def getTagPartnerTups(self, ldSNVtup):
        return self.ldPartnerTup2tagTups[ldSNVtup]

    def getTagSNVtups(self):
        return self.tagTup2ldPartnerTups.keys()


class GWAS_Annotations:
    """Contains annotation information for a single GWAS's tag SNPs and
    their LD partners."""

    def __init__(self, ldInfo, annotGroup):
        self.ldInfo = ldInfo
        self.annotGroup = annotGroup

    def writeToCSV(self, outfile, dnaseTissue = None):
        # For each tag SNV...
        for tagSNVtup in self.ldInfo.getTagSNVtups():
            # First column is the SNV's name and coords:
            outStr = tagSNVtup[0] + "_" + str(tagSNVtup[1]) + "_" + tagSNVtup[2]

            # FIXME: This still seems ugly, in terms of the representation
            # of linkage between SNPs. Will probably have to refactor
            # yet again at some point.

            # Get all LD SNP tups, in this form: ((chrom, loc), r2).
            ldSNVtups = self.ldInfo.getLDpartnerTups(tagSNVtup)

            # Adds a comma-separated list of their coords to the output string:
            ldSNVstrs = map(lambda tup: tup[0][0] + "_" + str(tup[0][1]) + str(tup[1]),
                            ldSNVtups)
            ldSNVtupStr = reduce(lambda tok1, tok2: tok1 + "," + tok2,
                                 ldSNVstrs)

            ldSNVlocs = map(lambda tup: (tup[0][0], tup[0][1]), ldSNVtups)
            currSNVannots = map(lambda ldSNVloc: self.annotGroup.getAnnotFromLoc(ldSNVloc),
                                ldSNVlocs)
            outStr = outStr + " " + combineAnnotsToCSV(currSNVannots, dnaseTissue = dnaseTissue)

            print >> outfile, outStr


def combineAnnotsToCSV(snvAnnots, dnaseTissue = None):
    """Computes some summary overlap stats over all the specified SNPs
    (e.g. "or" function on each of the features). Returns a string CSV
    representation."""

    # Currently:
    # Generate comma-separated lists of *number* of DNase overlap,
    # footprint overlap, and ChIP-seq nearby:

    # If told to look for a particular tissue when summarising DNase,
    # then do so:
    if (dnaseTissue == None):
        nDNaseArr = map(lambda snvAnnot: str(len(snvAnnot.getDNase())),
                        snvAnnots)
    else:
        nDNaseArr = map(lambda snvAnnot: str(int(dnaseTissue in snvAnnot.getDNase())),
                        snvAnnots)
    if len(nDNaseArr) > 0:
        dnaseStr = reduce(lambda tok1, tok2: tok1 + "," + tok2,
                          nDNaseArr)
    else:
        dnaseStr = ""

    nFootprintArr = map(lambda snvAnnot: str(len(snvAnnot.getFootprint())),
                        snvAnnots)
    if len(nFootprintArr) > 0:
        footprintStr = reduce(lambda tok1, tok2: tok1 + "," + tok2,
                              nFootprintArr)
    else:
        footprintStr = ""

    nChipseqArr = map(lambda snvAnnot: str(len(snvAnnot.getChipseq())),
                      snvAnnots)
    if len(nChipseqArr) > 0:
        chipseqStr = reduce(lambda tok1, tok2: tok1 + "," + tok2,
                            nChipseqArr)
    else:
        chipseqStr = ""

    # Just get the single *best* motif effect for each ld SNP. Generate
    # comma-separated list of those best ones:
    motifEffectArr = map(lambda snvAnnot: snvAnnot.getBestMotifEffect(),
                         snvAnnots)
    motifEffectStrs = map(lambda motEffect: outputHelperFn(motEffect),
                          motifEffectArr)
    if len(motifEffectStrs) > 0:
        motifEffectStr = reduce(lambda tok1, tok2: tok1 + "," + tok2,
                                motifEffectStrs)
    else:
        motifEffectStr = ""

    # Combine the lists of bed overlaps for the snvs:
    bedOverlapsArrs = map(lambda snvAnnot: snvAnnot.getBEDannots(),
                          snvAnnots)
    bedOverlapsCombined = reduce(lambda arr1, arr2: arr1+arr2, bedOverlapsArrs)
    if (len(bedOverlapsCombined) > 0):
        bedOverlapsStr = reduce(lambda str1, str2: str1 + "," + str2,
                                bedOverlapsCombined)
    else:
        bedOverlapsStr = ""

    # Report the *maximum* feature score and percentile, maximising over all ld
    # proxies, as a separate tuple for each feature name...
    # This was initially introduced to deal with DNase signal only (and may stay
    # that way).
    featureScoresAndPercentileTups = []
    for featureName in snvAnnots[0].name2scoreAndPercentile.keys():
        highestPercentile = -1
        highestScore = -1
        for snvAnnot in snvAnnots:
            try:
                if snvAnnot.name2scoreAndPercentile[featureName][1] > highestPercentile:
                    highestScore = snvAnnot.name2scoreAndPercentile[featureName][0]
                    highestPercentile = snvAnnot.name2scoreAndPercentile[featureName][1]
            except KeyError, e:
                # HACK: Dealing with weirdo SNPs; just ignoring for the time being
                # (September 25th 2014):
                pass

        currTupString = "%s/%d/%1.2f" % (featureName, highestScore, highestPercentile)
        featureScoresAndPercentileTups.append(currTupString)
    featureScoresAndPercentilesStr = "_".join(featureScoresAndPercentileTups)

    # DNase PWM:
    # Just report the maximum preferred allele minus log10 p-value, maximising
    # over all SNP/TF/PWM combinations:
    # FIXME: Not bothering to implement this yet.
    dnasePWMstr = ""

    # Combine into the resulting output string:
    return str(dnaseStr) + " " + str(footprintStr) + " " + str(chipseqStr) + " " + \
        str(motifEffectStr) + " " + bedOverlapsStr + " " + featureScoresAndPercentilesStr + \
        " " + dnasePWMstr


def outputHelperFn(motEffect):
    if motEffect == None:
        return "-"
    else:
        return motEffect.toString(outFormat="underscore")
