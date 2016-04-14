#!/usr/bin/env python

# A module for facilitating reading and writing of binary wiggle file data.


import commands, math, os, pdb, random, string, struct, sys


# Simple utility function for reversing elements of an input array:
def reverse(inArray):
    outArr = []
    idx = 1
    while (idx <= len(inArray)):
        outArr.append(inArray[-idx])
        idx = idx + 1
    return outArr


class wibConverter:
    """A class for converting between binary and floating representation of
    wiggle data."""

    def __init__(self, lower, upper, null):
        self.lower = lower
        self.upper = upper
        self.null = null

    def getRange(self):
        return self.upper-self.lower

    def binary2float(self, byte):
        # Convert binary byte value into float:
        if (byte == '\xFF'):
            return self.null
        else:
            # Calculate the ratio of the current value relative to the full
            # range; i.e. a value between 0 and 1 inclusive...
            byteIntVal = struct.unpack('B', byte)[0]
            rangeRatio = byteIntVal/float(254) # 254 comes from (byte - 2)

            # Calculate the wiggle value that this maps to, based on the
            # lower and upper values of this converter...
            fullRange = self.getRange()
            rangeVal = rangeRatio*fullRange + self.lower
            return rangeVal

    def float2binary(self, wigFloat):
        # Convert float into binary byte value:
        if (wigFloat == self.null):
            return '\xFF'

        if (self.lower > wigFloat):
            print >> sys.stderr, "WARNING: Value ", wigFloat, "was lower than",\
                "the specified lower bound value. Converting to", self.lower,\
                "instead."
            wigFloat = self.lower
        if (self.upper < wigFloat):
            print >> sys.stderr, "WARNING: Value ", wigFloat, "was higher",\
                "than the specified upper bound value. Converting to", \
                self.upper, "instead."
            wigFloat = self.upper

        # Calculate the fraction of the range that the wiggle value
        # represents:
        fullRange = self.getRange()
        rangeRatio = (wigFloat - self.lower)/float(self.getRange())
        
        # Calculate the byte value that this ratio maps to. This is where
        # the precision is lost:
        byteValFloat = rangeRatio*254
        byteValInt = int(round(byteValFloat)) # Ceil or floor instead??
        return struct.pack('B', byteValInt)


class wibWriter:
    """A writer for wib files."""

    def __init__(self, outDir, infile, conv):
        # Create the output directory if it does not already exist:
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        self.outDir = outDir
        self.wigInfile = infile
        self.conv = conv

        self.currPos = 1 # Curr position in chromosome (of elem being written).

        # Read in the current line from wigInfile and set self.currLine to it:
        self.currLine = infile.readline()

        # Parse the first line, setting info on the chromosome and curr block...
        elems = self.currLine.strip().split()

        # If the first line is a track description line then just ignore it:
        if (elems[0] == "track"):
            self.currLine = infile.readline()
            elems = self.currLine.strip().split()

        assert (elems[0] == "fixedStep")

        self.currBlockChrom = elems[1].split("=")[1]
        self.currBlockStart = int(elems[2].split("=")[1])
        self.currBlockStepSize = int(elems[3].split("=")[1])

        self.wibOutfile = open(self.outDir + '/' + self.currBlockChrom + '.wib', 'w')

    def writeWib(self):
        """Parses in the input wiggle file and writes corresponding output to
        the output binary wib file."""

        while (self.currLine != ""):
            # Process the current fixedStep block information:
            self.processBlock()
        self.wibOutfile.flush()

    def processBlock(self):
        # Write out null values preceding the current block...
        while (self.currPos < self.currBlockStart):
            self.wibOutfile.write(self.conv.float2binary(self.conv.null))
            self.currPos = self.currPos + 1

        # Skip over the block information line:
        self.currLine = self.wigInfile.readline()

        while ((self.currLine != "") and
               (self.currLine.split()[0] != "fixedStep")):
            # The current line designates a wiggle value for the current block.
            # Parse it and record the value...
            currVal = float(self.currLine)

            currValByte = self.conv.float2binary(currVal)

            # Write out currBlockStepSize instances of the byte to the binary
            # output file...
            for i in range(self.currBlockStepSize):
                self.wibOutfile.write(currValByte)
                self.currPos = self.currPos + 1
            self.currLine = self.wigInfile.readline()

        # Parse and store the block information for the next block:
        if (self.currLine != ""):
            elems = self.currLine.strip().split()

            assert (elems[0] == "fixedStep")

            self.currBlockChrom = elems[1].split("=")[1]
            self.currBlockStart = int(elems[2].split("=")[1])
            self.currBlockStepSize = int(elems[3].split("=")[1])


class wibReader:
    """A reader for wib files."""

    def __init__(self, inDir):
        self.inDir = inDir
        rangeFile = open(inDir + "/range.txt")
        lower = float(rangeFile.readline().strip())
        upper = float(rangeFile.readline().strip())
        null = float(rangeFile.readline().strip())
        conv = wibConverter(lower, upper, null)
        self.upper = upper
        self.lower = lower
        self.null = null
        self.conv = conv

    def getNull(self):
        return self.null

    def getWigs(self, track, flipOnStrand=False):
        """Retrieves the wiggle data for the specified BED_Track, returning
        it as a dictionary where each key is a BED_Line and each value is
        an array of the wiggle data."""

        wigData = {}
        for bedLine in track.getAllBEDs():
            wigData[bedLine] = self.getWig(bedLine, flipOnStrand)

        return wigData

    def getWig(self, line, flipOnStrand=False):
        """Retrieves wiggle data for the single specified line, and records
        that information internally."""

        # Extract the chromosome, start, and end for the region...
        chrom = line.getChrom()
        start = line.getChromStart()
        end = line.getChromEnd()

        # Check that there is a file for that chromosome in the specified
        # inputDir; report a warning and return if there is not...
        wibFilename = self.inDir + '/' + chrom + '.wib'
        if (not os.path.exists(wibFilename)):
            print >> sys.stderr, "WARNING: No chromosome file found for " + \
                chrom
            return

        # Open the chromFile:
        chromWibFile = open(wibFilename)

        # Skip over bytes in the file until the specified start is
        # reached:
        try:
            chromWibFile.seek(start-1)

            # Read the data for the current block of interest (end-start bytes):
            wibData = chromWibFile.read(end-start+1)

            # Convert the data to wiggle format, from binary wib format...
            wigData = []

            for byte in wibData:
                wigData.append(self.conv.binary2float(byte))

            # Consider strandedness:
            if flipOnStrand:
                if (line.getStrand() == "+"):
                    return wigData
                else:
                    assert (line.getStrand() == "-")
                    return reverse(wigData)
            else:
                return wigData
        except Exception, e:
            print >> sys.stderr, start, chromWibFile
#            raise e
