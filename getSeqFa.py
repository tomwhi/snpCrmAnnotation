#!/usr/bin/env python

import commands, math, os, pdb, re, string, sys
import genomics.formats.bed as bed
import seqUtils
from pylab import *
from optparse import OptionParser

# Note: I'm currently assuming the newline character is 1 long:
newLineLen = 1

class seqRetriever:
    """This class is used in order to retrieve sequences of specified
    coordinates from a specified fasta file.
    """

    def __init__(self, faFilename):
        self.seqFilename = faFilename

        # Create a new index of the fasta file if necessary:
        self.indexFilename = faFilename + ".fai"
        # Check whether a corresponding index file exists:
        if (not (os.access(self.indexFilename, os.F_OK))):
            # There is no pre-existing index file => create one:
            indexer = faFileIndexer(self.seqFilename)
            indexer.extractIndex()
            indexFile = open(self.indexFilename, 'w')
            indexer.writeIndexToFile(indexFile)
            indexFile.flush()
            indexFile.close()

        
        # Parse the found (or created) index file into an index dictionary,
        # storing that index as an attribute of this object:
        newIndexer = faFileIndexer()
        newIndexer.parseIndexFile(open(self.indexFilename))
        self.index = newIndexer.getIndex()

    def getSequence(self, loc, revComp=False):
        """Extract the sequence corresponding to the specified location, from
        the sequence of this seqRetriever.

        Algorithm:
        - Retrieve the index tuple corresponding to specified sequence name,
          from the faIndex dictionary.
        - Calculate the start and end positions of the required sequence in the
          input file (i.e. number of characters preceding that position in the
          file), using the following formula:

          pos = r + x + h*q - 1
          q = floor((x-1)/p)

          where:
          pos is the number of characters preceding the required position in
          the file
          r is the number of characters preceding the desired sequence
          x is the position in the sequence, counting from 1
          p is the number of sequence characters on each line
          q is the number of overhang lines
          h is the number of overhang characters on each line (including any
          newline characters)

        - Perform a "seek" to the starting position in the input fasta file,
          and read the end-start characters into a string.
        - Split that string into an array of occurrences of any newline
          instances. Return this array."""

        # Check precondition that the location matches "<chr>:<start>-<end>"...
        if (re.search('.*:.*-.*', loc) == None):
            print >> sys.stderr, "Specified chromosomal location is invalid: ",\
                loc
            sys.exit(1)

        # Parse the specified location:
        toks1 = loc.split(":")
        if not (len(toks1) == 2):
            print >> sys.stderr, "Specified chromosomal location is invalid: ",\
                loc
            sys.exit(1)
        seqName = toks1[0]
        toks2 = toks1[1].split("-")
        if not (len(toks2) == 2):
            print >> sys.stderr, "Specified chromosomal location is invalid: ",\
                loc
            sys.exit(1)
        start = int(toks2[0])
        end = int(toks2[1])
        
        # Prepare to read from the sequence fasta file:
        faFile = open(self.seqFilename)

        # Retrieve the index tuple corresponding to the specified sequence name,
        # and the required values from that index:
        if (not (self.index.has_key(seqName))):
            # The inputs are inconsistent => Report this by raising an
            # exception:
            raise ValueError("No sequence with name " + seqName +
                             " in file " + self.seqFilename + ".")

        seqIndex = self.index[seqName]
        h = seqIndex[3] - seqIndex[2]
        p = seqIndex[2]
        r = seqIndex[1]

        # Calculate the starting position for the file read:
        fileStartPos = int(r + start + h*floor((start-1)/float(p)) - 1)
        
        # Perform a seek to that starting position:
        faFile.seek(fileStartPos)

        # Read end-start characters into a string:
        fileEndPos = int(r + end + h*floor((end-1)/float(p)) - 1)

        # Include both the start and end positions:
        nCharsToRead = fileEndPos - fileStartPos + 1
        outStr = faFile.read(nCharsToRead)

        # Remove any newline characters:
        outSeq = outStr.replace("\n","")

        # Return the resulting sequence, or its reverse complement if requested:
        if (revComp):
            return seqUtils.revcomp(outSeq)
        else:
            return outSeq


class faFileIndexer:
    """An object of this class is used to generate an index file for a single
    fasta file."""

    def __init__(self, infileName=None):
        self.faFile = None
        if (infileName != None):
            # This index is going to be created from a specified input file:
            self.faFile = open(infileName)

        self.nCharsRead = 0
        self.index = {}
        self.currLine = None

    def extractIndex(self):
        # Extract the index from the input file...

        # Must have not yet read anything from the input file:
        assert (self.nCharsRead == 0)
        assert (len(self.index) == 0)

        # Keep reading and indexing sequences from the fasta file until
        # the entire fasta file has been read...

        self.currLine = self.faFile.readline()
        while (self.currLine != ""):
            self.extractSeqIdx()
    
    def extractSeqIdx(self):
        """Parses a single sequence from a specified input stream, and adds
        an entry to the index, for that sequence."""

        currSeq_lineLen = None
        currSeq_overhang = None

        # First, read the sequence name line, calculate the length of the name
        # line, and extract the sequence name as a key for the index:
        assert (self.currLine[0] == ">")
        seqName = self.currLine[1:].strip()

        # Make sure that the sequence name adheres to the file format
        # requirements:
        assert (len(seqName) > 0)
        assert (not (self.index.has_key(seqName)))

        # Update the current position within the file.
        self.nCharsRead = self.nCharsRead + len(self.currLine)

        # That position is the position at which this sequence starts. Record
        # this, to be written to this sequence's index information:
        seqStartPos = self.nCharsRead

        nSeqLinesRead = 0
        self.currLine = self.faFile.readline()
        while ((self.currLine != "") and (self.currLine[0] != ">")):
            # Read a single fragment of sequence...

            if (currSeq_lineLen == None):
                # Determine the length of sequence lines in this sequence, in
                # order to add that info to this sequence's index:
                currSeq_lineLen = len(self.currLine) - newLineLen

                # Presently, setting the overhang amount statically, assuming
                # it to be "1" (for one newline character):
                currSeq_overhang = newLineLen

            # Otherwise, no more index information needs to be gleaned from
            # this line, although the position of the reader within the
            # file needs to be updated. NOTE: I could in principle
            # implement some error checking, e.g. making sure the input
            # lines are all the same length except for the last one. Presently
            # I am not doing this error checking.

            # Update this reader's current position within the fasta file:
            self.nCharsRead = self.nCharsRead + len(self.currLine)

            # Record the fact that a sequence line was read:
            nSeqLinesRead = nSeqLinesRead + 1

            self.currLine = self.faFile.readline()

        # Record the index information for the current sequence (which is
        # [seqLen, nPrecedingChars, seqLineLen, totalLineLen]
        seqLen = self.nCharsRead - (nSeqLinesRead*newLineLen) - seqStartPos
        self.index[seqName] = [seqLen, seqStartPos, currSeq_lineLen,
                               currSeq_lineLen + newLineLen]


    def writeIndexToFile(self, outfile):
        """Writes an index for <faFilename> to the specified output file,
        from the index generated via a call to this object's "extractIndex()"
        method."""

        for seqName in self.index.keys():
            outfile.write(seqName + "\t")
            indexTuple = self.index[seqName]
            outfile.write(reduce(lambda tok1, tok2: str(tok1) + "\t" + \
                                     str(tok2), indexTuple) + "\n")

    def parseIndexFile(self, indexFile):
        """Generate this indexer by reading an index from the specified file."""
        for line in indexFile:
            elems = line.strip().split("\t")
            seqName = elems[0]
            seqIndexToks = elems[1:]
            seqIndex = [int(seqIndexToks[0]), int(seqIndexToks[1]),
                        int(seqIndexToks[2]), int(seqIndexToks[3])]
            self.index[seqName] = seqIndex

    def getIndex(self):
        return self.index


def main():
    # Parse the command-line arguments...
    description = """usage: %prog [options] <faFile>\n

Retrieves sequences corresponding to the genomic regions specified as a bed file
on stdin. The genomic sequence file is specified as a command-line argument.
Uses indexing to do this quickly. If an index file (of the name <faFile>.fai)
exists already, it just uses that index. Otherwise it generates a new index
file. This program is similar to samtool's faidx utility.
"""

    parser = OptionParser(usage = description)
    parser.add_option("--debug", action="store_true", dest="debug",
                      help = "Debug the program using pdb.")
    parser.add_option("--considerStrand", dest = "considerStrand",
                      default = False, action = "store_true",
                      help = "Consider strand of each input bed item, and " + \
                          "print out reverse complement of the sequence " + \
                          "when region is on minus strand.")
    (options, args) = parser.parse_args()

    if (options.debug):
        pdb.set_trace()

    # Parse the input parameters...

    # Make sure the required input file exists and that the seqRegion
    # parameter was specified correctly:
    if (len(args) != 1):
        print >> sys.stderr, "WRONG # ARGS: ", len(args)
        parser.print_help()
        sys.exit(1)

    # Prepare to extract the specified sequence:
    seqGetter = seqRetriever(args[0])

    # Retrieve the requested sequences:
    for line in sys.stdin.xreadlines():
        lineAsBed = bed.BED_line(line.strip())
        regStr = lineAsBed.getChrom() + ":" + str(lineAsBed.getChromStart()) + \
            "-" + str(lineAsBed.getChromEnd())
        outSeq = None
        if (options.considerStrand):
            if (lineAsBed.getStrand() == "+"):
                outSeq = seqGetter.getSequence(regStr)
            else:
                assert (lineAsBed.getStrand() == "-")
                outSeq = seqGetter.getSequence(regStr, revComp=True)
        else:
            outSeq = seqGetter.getSequence(regStr)

        print ">" + lineAsBed.to_string()
        print outSeq

if __name__ == "__main__":
    sys.exit(main())



#<seqRegion> must follow the format:
#<seqName>:<start>-<end>
#, where <seqName> is the name of a sequence present within <faFile>, and the
#<start>/<end> coordinates (counting from [1]) are contained within the span
#of that sequence.
