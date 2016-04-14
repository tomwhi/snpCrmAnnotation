#!/usr/bin/env python

# A module for generating wiggle files (e.g. genome tag-count coverage scores
# from bam and bed files.

import gc, os, pdb, sys
import arrayUtils, cmdlineProgs, utility
import xml.etree.ElementTree as ET
import genomics.formats.bed as bed
import genomics.formats.conversion as conversion
import genomics.formats.wibIO as wibIO
#from guppy import hpy


def dict2varStep(chrom2pos2score, wigOutputFile, trackName="trackName",
                 trackDesc = "Track description"):
    """Utility function, generates a variable-step wiggle track
    from the specified chrom -> pos -> score dictionary."""

    # Print a track definition line:
    print >> wigOutputFile, "track name=\"" + trackName + "\" type=wig"

    # For each chromosome...
    for chrom in chrom2pos2score.keys():
        # Write a block header for this chromosome:
        print >> wigOutputFile, "variableStep chrom=" + chrom + " span=1"

        # Will report the scores in order of chromosomal coordinates:
        currChromPositions = chrom2pos2score[chrom].keys()
        currChromPositions.sort()
        
        # For each position for this chromosome...
        for pos in currChromPositions:
            # Print the coordinate and the score:
            print >> wigOutputFile, pos, chrom2pos2score[chrom][pos]


class Profile:
    """Represents a "wiggle" track profile and the set of files (wig/wib/tdf)
    storing that profile."""

    def __init__(self, name, prefix, params, genome):
        # A name for the profile, used when generating header lines.
        self.name = name

        # The start of location of the files representing this profile:
        # FIXME: Change this way of doing things if I want to make the
        # tdf/wig/wib file locations more flexible.
        self.prefix = prefix

        # A ProfileParams object defining the profile (e.g. extension length,
        # combined/subtracted/plus-strand/minus-strand):
        self.params = params

        self.wigExists = False
        self.wibExists = False
        self.tdfExists = False

        self.genome = genome

    def getWigLoc(self):
        return self.prefix + ".wig"

    def getWibLoc(self):
        return self.prefix + "_wib"

    def getTDFLoc(self):
        return self.prefix + ".tdf"

    def getHeaderName(self):
        """Returns a string representation to act as header names for this profile"""
        return self.name + "_" + self.params.toString()

    def deleteWig(self):
        """Deletes the wiggle file for this profile."""

        cmdlineProgs.deleteFiles([self.getWigLoc()])
        self.wigExists = False

    def makeTDF(self):
        """Converts the wig file from wiggle to tdf format, keeping the original
        file."""

        conversion.wig2tdf(self.getWigLoc(), self.getTDFLoc(),
                           genome=self.genome)
        self.tdfExists = True

    def makeWib(self, lower = 0.0, upper = 1.0, null = -1.0):
        """Converts the wig file from wiggle to wib format, keeping the original
        file."""

        # FIXME: Currently, this produces a directory, not a single wib file.
        # Potential Solution: Should I use the bigWig format instead, and just
        # stop using my own wib format?? -> I won't yet, I have a quite nice
        # interface to the wib format (e.g. it is already integrated with the
        # bed.BED_Track class).
        self.wibExists = True

        conv = wibIO.wibConverter(lower, upper, null)

        # Set up a writer to write the wib data:
        wibWriter = wibIO.wibWriter(self.getWibLoc(), open(self.getWigLoc()), conv)

        # Write the binary wib file:
        wibWriter.writeWib()

        # Write out the range file specifying the lower, upper, and null values:
        rangeFile = open(self.getWibLoc() + "/range.txt", 'w')
        print >> rangeFile, lower
        print >> rangeFile, upper
        print >> rangeFile, null

    def editXMLelem(self, xmlElem):
        """Adds an XML summary of this profile to the specified XML document."""

        # Store the locations of the wib and tdf files:
        profileNode = ET.SubElement(xmlElem, 'Profile')
        profileNode.text = self.getHeaderName()
        if self.wibExists:
            wibDirNode = ET.SubElement(profileNode, 'WibDir')
            wibDirNode.text = self.getWibLoc()
        if self.tdfExists:
            tdfFileNode = ET.SubElement(profileNode, 'TDFFile')
            tdfFileNode.text = self.getTDFLoc()

    # July 25th 2014: Stopped using this.
    def writeIGV_XML(self, parentNode, baseDataDir2urlBase, nameForXML, att2val={}):
        """Generate an XML element given a parent XML element. This is
        for the IGV server for displaying the TDF track of this profile."""

        # Precondition:
        assert self.tdfExists

        currParamsNode = ET.SubElement(parentNode, "Resource")
        currParamsNode.set("name", nameForXML)
        for att in att2val:
            currParamsNode.set(att, att2val[att])

        # Determine TDF URL:
        tdfLocAbs = os.path.abspath(self.getTDFLoc())
        urlLoc = tdfLocAbs.replace(baseDataDir2urlBase[0], baseDataDir2urlBase[1])

        # Record the TDF URL location in the XML:
        currParamsNode.set("path", urlLoc)

    # FIXME: Later, implement the option to specify the resolution of the
    # wiggle data (by performing some averaging/windowing procedure); currently
    # it is fixed at 1bp.


class ProfileParams:
    """Represents the parameters for a single wiggle profile: Extension length,
    and type (one of "plus", "minus", "combined", "subtracted")."""

    plus = "plus"
    minus = "minus"
    combined = "combined"
    subtracted = "subtracted"

    def __init__(self, extension, profType):
        self.extension = extension
        self.profType = profType

    def toString(self):
        return self.profType + "_" + str(self.extension)

    def getExtension(self):
        return self.extension

    def getProfType(self):
        return self.profType


def getBlockCoords(bedTrack, maxRegs=50):
    """Returns a list of chromosomal coordinate strings for blocks encompassing
    the specified set of bed items, with up to 500 bed items per block.

    This is necessary in order to make profile generation fast and also
    not consume too much memory."""

    # The list of coordinates strings to be returned:
    blockCoords = []

    # Keep each block within a single chromosome...
    for chrom in bedTrack.get_chroms():
        # Get all regions for that chromosome:
        currChromBEDs = bedTrack.get_BEDs(chrom)

        # Split those regions into blocks...
        currChromRegsBlocks = arrayUtils.divideList(currChromBEDs, maxRegs)


        # For each of those blocks...
        for currBlockRegs in currChromRegsBlocks:
            # Get the chrom, start of the first element and the chrom, end of
            # the last element:
            firstReg = currBlockRegs[0]
            lastReg = currBlockRegs[-1]

            blockStart = firstReg.getChromStart()
            blockEnd = lastReg.getChromStart()

            # Generate a coordinate string from that chrom, start, end:
            currBlockCoords = chrom + ":" + str(blockStart) + "-" + str(blockEnd)

            blockCoords.append(currBlockCoords)

    return blockCoords


class ProfileGenerator:
    """Used for generating a wiggle profile given an input bam and bed file and
    parameters specifying how to generate the profile."""

    def __init__(self, name, prefix, paramsForProfs, genome,
                 minScore=-127, maxScore=128):
        self.name = name
        self.prefix = prefix
        self.genome = genome
        self.minScore = minScore
        self.maxScore = maxScore

        # For retrieving params objects given profile type and extension:
        self.paramsTup2params = {}

        # Dictionary of profileParams->profile objects. Stores the wiggle
        # profiles and their parameters, as generated by this profile generator:
        self.params2prof = {}
        for params in paramsForProfs:
            currProfPrefix = self.prefix + "_" + params.toString()
            self.params2prof[params] = \
                Profile(self.name, currProfPrefix, params, self.genome)

            self.paramsTup2params[(params.getProfType(), params.getExtension())] = \
                params

    def getAllParams(self):
        return self.params2prof.keys()

    def getAllProfs(self):
        return self.params2prof.values()

    def getProf(self, profType, extension):
        params = self.paramsTup2params[(profType, extension)]
        return self.params2prof[params]

    def makeAllRegProfiles(self, bamFilename, regsBedTrack, verbose=False):
        """Generates all wiggle profiles for each of the parameter sets for
        this generator. Writes the .wig files, converts to .tdf and .wib and
        deletes the .wig files."""

        # If there are no profiles specified by the parameters, then
        # don't do anything here:
        if len(self.params2prof.values()) == 0:
            print >> sys.stderr, \
                "Not generating profiles since no profile parameters were specified."
            return

        # Start out all wig files for the profiles, and write their
        # headers:
        # FIXME: Generating the wig files in parts by opening them each
        # time. Must make sure to flush the files after each write, otherwise
        # the file could be jumbled. It might be more elegant to keep a
        # reference to the open (writing) wig file, in the profile object.
        for profile in self.params2prof.values():
            wigFile = open(profile.getWigLoc(), 'w')
            headerLine = "track name=" + profile.getHeaderName() + \
                " description=\"" + profile.getHeaderName() + "\""
            print >> wigFile, headerLine
            wigFile.close()

        # Temporary file for containing reads for current region, in
        # bed format:
        readsTmpFilenameBed = utility.makeTempFilename("CurrChrom_bed")

        # Reads with name set to X (to save space when loaded to memory):
        readsNoName_Filename = utility.makeTempFilename("CurrChrom_bed_NoName")

        # For each currRegion in the bed file, write out a profile for it...

        # August 1st 2013: Process the regions in blocks of (default) 50,
        # to balance memory usage and time usage:
        blockCoords = getBlockCoords(regsBedTrack, maxRegs=50)

        # Process "profile elements" (which will be used to build the profiles
        # e.g. sequence reads) one block at a time:
        for currBlockCoords in blockCoords:
            # The block bed tracks can be gigantic; need to collect unused
            # ones to make sure memory usage doesn't blow up:
            currBlockBEDs = None
            currRegionBEDs = None
            unreachable=gc.collect()
            print >> sys.stderr, "Ran gc.collect(). Result:", unreachable

            #h = hpy()
            #print >> sys.stderr, "Heap info from guppy:"
            #print >> sys.stderr, h.heap()

            # Truncate negative coords to 1:
            blockCoords = currBlockCoords.split(":")[1]
            if blockCoords[0] == "-":
                newBlockCoords = "1-" + blockCoords.split("-")[-1]
                currBlockCoords = currBlockCoords.split(":")[0] + ":" + newBlockCoords

            # Get all reads for the current chromosome, as a bed track:
            conversion.bam2bed(bamFilename, readsTmpFilenameBed,
                               regionStr=currBlockCoords, verbose=verbose)

            readsNoName_file = open(readsNoName_Filename, 'w')

            for line in open(readsTmpFilenameBed).readlines():
                elems = line.strip().split()
                print >> readsNoName_file, \
                    elems[0], elems[1], elems[2], "X", elems[4], elems[5]
            
            currBlockBEDs = bed.BED_Track(open(readsTmpFilenameBed))

            print >> sys.stderr, "Generating wiggle profiles for block", \
                currBlockCoords, "..."
            profIdx = 0

            # Process all regions in the current block:
            currBlockAsBEDline = \
                bed.BED_line(currBlockCoords.replace(":", " ").replace("-", " "))

            bedsForRegion = None
            try:
                bedsForRegion = regsBedTrack.getBEDsForRegion(currBlockAsBEDline)
            except Exception, e:
                pdb.set_trace()
                dummy = 1
            for currRegion in bedsForRegion:
                # Get all beds to appear in the profile for the current region:
                try:
                    currRegionBEDs = currBlockBEDs.getBEDsForRegion(currRegion)
                except Exception, e:
                    pdb.set_trace()
                    dummy = 1

                # Make profiles for the region:
                self.makeRegProfiles(currRegion, currRegionBEDs,
                                     verbose=verbose)
                if (profIdx % 10) == 0:
                    print >> sys.stderr, "  Progress:", str(profIdx), "regions."
                profIdx += 1
            print >> sys.stderr, "Done."

        #cmdlineProgs.deleteFiles([readsTmpFilenameBed, currChromBamFilename])

        # For each profile for this profile generator, convert to tdf and wib:
        for profile in self.params2prof.values():
            print >> sys.stderr, "Generating TDF for current profile..."
            profile.makeTDF()
            print >> sys.stderr, "Done."

            # FIXME: Currently, not generating a wib directory, since it would
            # consume way too much space (unless I implement a proper
            # fixedStep-indexed version inside wibIO.py).

            #print >> sys.stderr, "Generating wib directory for current profile..."
            #profile.makeWib(lower=self.minScore, upper=self.maxScore,
            #                null=self.minScore - 1)
            #print >> sys.stderr, "Done."

            # Delete the (potentially large, non-binary) wig file too:
            #profile.deleteWig()

    def makeRegProfiles(self, currRegion, bedsForRegion, verbose=True):
        """Generates all profiles for this profile generator, for a single
        region, and concatenates the info to the respective profile's wig
        files. bedsForRegion is a list of bed items constituting the
        profile elements, which reside in the current region."""

        if verbose:
            print >> sys.stderr, "Generating wiggle profile data for region" + \
                " " + currRegion.getBedCoordStr()

        # Get the reads for the current region, in bed format:
        regionStr = currRegion.getBedCoordStr()

        # For each of the parameter combinations...
        for params in self.params2prof:
            currProfile = self.params2prof[params]

            # Generate an array of the current region's profile given those
            # parameters and the bed file of reads for the region:
            currProfArr = self.profArrFromBEDs(params, currRegion, bedsForRegion)

            # Update the wiggle file of the corresponding profile with
            # that block of coverage scores:
            self.writeToWig(currProfArr, currRegion, currProfile.getWigLoc())

        if verbose:
            print >> sys.stderr, "Done."

    def profArrFromBEDs(self, params, region, bedsForRegion):
        """Generates an array of overlap scores for the specified region, given
        the profile parameters."""

        # region is a BED_line object.

        # Generate empty "scores" array for the region, starting as a zeros
        # array:
        scoresArr = [0] * region.getLen()

        # Parse the specified beds file, to get all bed item start positions
        # and strands, and increment/decrement counts in the two arrays
        # appropriately for each such bed item, taking into consideration the
        # extension length and "type" of the specified params...
        for bedItem in bedsForRegion:
            # Genomic position of the strand-specific start of the bed item:
            bedStart = None
            if (bedItem.getStrand() == "+"):
                bedStart = bedItem.getChromStart()
            else:
                bedStart = bedItem.getChromEnd()

            # Find the start and end coordinates (inclusive) of the (region of
            # the scores array that needs to be incremented) for the current bed
            # item:
            scoresStartPos = None
            if (bedItem.getStrand() == "+"):
                scoresStartPos = bedStart - region.getChromStart()
                scoresEndPos = scoresStartPos + params.extension - 1
            else:
                scoresEndPos = bedStart - region.getChromStart()
                scoresStartPos = scoresEndPos - params.extension + 1

            # Trim to ensure the region is inside the scores array:
            if scoresStartPos < 0:
                scoresStartPos = 0
            if scoresEndPos >= len(scoresArr):
                scoresEndPos = len(scoresArr) - 1

            # Determine the value to be added (+1, 0, or -1):
            toAdd = None
            if (params.profType == ProfileParams.plus):
                if (bedItem.getStrand() == "+"):
                    toAdd = 1
                else:
                    toAdd = 0
            elif (params.profType == ProfileParams.minus):
                if (bedItem.getStrand() == "-"):
                    toAdd = -1
                else:
                    toAdd = 0
            elif (params.profType == ProfileParams.combined):
                toAdd = 1
            else:
                assert (params.profType == ProfileParams.subtracted)
                if (bedItem.getStrand() == "+"):
                    toAdd = 1
                else:
                    toAdd = -1

            # Update the scores array at the calculated positions:
            scoresIdx = scoresStartPos
            while scoresIdx <= scoresEndPos:
                scoresArr[scoresIdx] = scoresArr[scoresIdx] + toAdd
                scoresIdx = scoresIdx + 1

        return scoresArr

    def writeToWig(self, profArr, region, wigFilename,
                   wigFormat="fixedStep"):
        """Writes the specified wiggle file segment out to the specified
        file."""
            
        # Header line:
        wigFile = open(wigFilename, 'a')

        # Currently only writing in fixedStep format:
        assert (wigFormat == "fixedStep")
        self.writeArrToFixedStep(profArr, region, wigFile)

        wigFile.close()

    def writeArrToFixedStep(self, profArr, region, wigFile):
        """Writes out the block of scores as a series of
        fixedStep wiggle format blocks."""

        # Break the scores array into a list of tuples, each tuple
        # specifying the score and the length it covers:
        lenAndScoreTups = arrayUtils.values2lenAndValueTups(profArr)

        # Merge each block of (consecutive tuples of the same length) into a
        # single tuple specifying that length and the array of tuples in
        # that block:
        lenAndScoreArrTups = \
            arrayUtils.tagAndValueTups2tagAndValueArrTups(lenAndScoreTups)

        # Write each of those (length, [scoreArray]) blocks out, each one
        # as a "fixedStep" block:
        currRegStart = region.getChromStart()
        for lenScoreArrTup in lenAndScoreArrTups:
            length = lenScoreArrTup[0]
            scores = lenScoreArrTup[1]
            currRegionHdr = "fixedStep chrom=" + region.getChrom() + \
                " start=" + str(currRegStart) + " step=" + str(length) + \
                " span=" + str(length)
            print >> wigFile, currRegionHdr
            for score in scores:
                print >> wigFile, score
            currBlockTotalLen = length*len(scores)
            currRegStart += currBlockTotalLen


def bigWig2percentiles(bigWigFilename, assumeFromZeroUnmeasured=True):
    """Processes a big wig file and returns a dictionary mapping from value
    to a percentile value.

    assumeFromZeroUnmeasured is true if zero scores exist but are
    not reported in the bigWig file - such as in the ENCODE DNase data."""

    # First, convert to a temporary wiggle file:
    tmpWigFilename = utility.makeTempFilename("TemporaryWiggleFile")
    cmdlineProgs.runBigWigToWig(bigWigFilename, tmpWigFilename, verbose=True)

    # Start out a dictionary of all values observed and their counts,
    # score2count:
    score2count = {}

    # Parse the wiggle file to extract score counts and to get percentiles...
    tmpWigFile = open(tmpWigFilename)
    nextLine = tmpWigFile.readline()
    while nextLine != "":
        (currValsObserved, nextLine) = processWiggleBlock(tmpWigFile, nextLine)
        for val in currValsObserved.keys():
            # Update scoreCounts with those values:
            if score2count.has_key(val):
                score2count[val] += currValsObserved[val]
            else:
                score2count[val] = currValsObserved[val]

    cmdlineProgs.deleteFiles([tmpWigFilename])

    # Generate the score->percentile mapping...
    score2percentile = {}
    scoresSorted = score2count.keys()
    scoresSorted.sort()
    totalCount = reduce(lambda count1, count2: count1 + count2, score2count.values())
    if assumeFromZeroUnmeasured:
        # Need to correct for unreported zero scores, by assuming that these are
        # the lowest scores and occur at every position in the genome not
        # encompassed by reported scores...

        # Calculate the zero-score percentiles, by considering the total size of the genome:
        # FIXME: HACK: Hard-coding hg19 genome size here:
        totalGenomeSize = 3137144693 # From http://www.ncbi.nlm.nih.gov/assembly/2758/

        numZeroScores = totalGenomeSize - totalCount
        score2percentile[0] = float(numZeroScores)/totalGenomeSize

        currCountSum = numZeroScores
        for score in scoresSorted:
            currCountSum += score2count[score]
            score2percentile[score] = float(currCountSum)/totalGenomeSize
    else:
        currCountSum = 0
        for score in scoresSorted:
            currCountSum += score2count[score]
            score2percentile[score] = float(currCountSum)/totalCount

    return score2percentile


def processWiggleBlock(wigFile, nextLine):
    """Parses a single "variableStep" block, recording a count of how many times
    each value is observed. Returns a dictionary with those value->count
    pairs."""

    # Determine the currSpan value by parsing nextLine:
    assert nextLine[:12] == "variableStep"
    elems = nextLine.strip().split()
    span = int(elems[-1].split("=")[1])

    val2count = {}
    nextLine = wigFile.readline()

    #print >> sys.stderr, "PROCESSING WIGGLE BLOCK:", nextLine
    #print >> sys.stderr, "FIRST ELEMENT:", nextLine
    while nextLine != "" and nextLine[:12] != "variableStep":
        elems = nextLine.strip().split()
        currValue = int(elems[1])
        if val2count.has_key(currValue):
            val2count[currValue] += span
        else:
            val2count[currValue] = span
        nextLine = wigFile.readline()

    return (val2count, nextLine)


def parseBigWigAvgOutput(bigWigAvgFile):
    """Parses a tab-delimited output file from a bigWigAverageOverBed
    command, and returns a snpLoc2score dictionary."""
    
    snpLoc2bigWigScore = {}
    
    for line in bigWigAvgFile.xreadlines():
        elems = line.strip().split("\t")
        
        # Parse the snp location:
        snpChrom = elems[0]
        snpPos = int(elems[1])
        
        # FIXME: Currently assuming discrete integer scores, to facilitate
        # percentile value retrieval easily:
        meanScore = int(elems[-1])
        
        # Update the dictionary with this pair:
        snpLoc2bigWigScore[(snpChrom, snpPos)] = meanScore
        
    return snpLoc2bigWigScore
