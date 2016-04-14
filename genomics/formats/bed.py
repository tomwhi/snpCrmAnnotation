#!/usr/bin/env python

import sys
import math, pdb, random
import list_funcs
import array
#import pysam


class BED_Track:
    """A BED track, following the format defined in the UCSC browser website."""
    def __init__(self, trackInput, peaksTrack=False, refGenomeFilename=None):
        # Each chromosome of BED items will be held in an array, referenced
        # by a dictionary. This allows separate chromosome's data to be
        # easily retrieved:
        self.chroms = {}

        # Determine the type of input, and initialise this BED_Track according
        # to that type:
        if (isinstance(trackInput, file)):
            self.initFromFile(trackInput, peaksTrack=peaksTrack)
        else:
            assert (isinstance(trackInput, list))
            self.initFromList(trackInput)

        # The name of a file containing the reference genome corresponding
        # to the assembly for this BED_Track object.
        self.refGenomeFilename = refGenomeFilename

    def getRefGenomeFilename(self):
        return self.refGenomeFilename

    def setRefGenomeFilename(self, filename):
        self.refGenomeFilename = filename

    def initFromList(self, inList):
        """Private method for initialising this object from an input list."""
        for inputBED in inList:
            assert(isinstance(inputBED, BED_line))
            bed_chrom = inputBED.getChrom()
            # Add the bed data item to the correct chromosome:
            if self.chroms.has_key(bed_chrom):
                self.chroms[bed_chrom].append(inputBED)
            else:
                self.chroms[bed_chrom] = [inputBED]

        # Sort the bed_items according to chromosomal location for each
        # chromosome:
        for key in self.chroms.keys():
            bed_items = self.chroms[key]
            bed_items.sort(chrom_compare)
            tuple_bed_items = tuple(bed_items)
            self.chroms[key] = tuple_bed_items

    def initFromFile(self, inFile, peaksTrack=False):
        """Private method for initialising this object from an input file."""
        # Parse the BED track from the specified file, creating a BED_Track
        # object:
        try:
            for line in inFile:
                elems = line.split()
                # Only add the line if it is a chromosomal location:
                #if ((len(elems) > 0) and (elems[0][0:3] == "chr")):
                # FIXME: The above condition was broken, so I changed it
                # to only filter out comment lines:
                if ((len(elems) > 0) and (elems[0][0] != "#")):
                    # Generate either a regular BED_line or a peaks_BED_line,
                    # depending on the value of the optional peaksTrack
                    # parameter. This is done to allow this class to optionally
                    # represent a track of ChIP-seq peak data, with extra
                    # attributes on each line:
                    curr_bed = None
                    if (peaksTrack):
                        curr_bed = peak_BED_line(line, track=self)
                    else:
                        curr_bed = BED_line(line, track=self)

                    bed_chrom = curr_bed.getChrom()
                    # Add the bed data item to the correct chromosome:
                    if self.chroms.has_key(bed_chrom):
                        self.chroms[bed_chrom].append(curr_bed)
                    else:
                        self.chroms[bed_chrom] = [curr_bed]
            # Sort the bed_items according to chromosomal location, if
            # necessary, for each chromosome:
            for key in self.chroms.keys():
                bed_items = self.chroms[key]
                bed_items.sort(chrom_compare)
                tuple_bed_items = tuple(bed_items)
                self.chroms[key] = tuple_bed_items
                
        except IOError:
            # If the file does not open, throw this exception upwards, rather
            # than dealing with it here:
            raise

        except ValueError:
            # A value error here means that a BED_line was invalid.
            pdb.set_trace()
            dummy = 1
            raise

    # Added on 10th November 2009:
    def getAllBEDs(self):
        allBEDs = []
        for chrom in (self.get_chroms()):
            allBEDs = allBEDs + list(self.get_BEDs(chrom))
        return allBEDs

    # Function to retrieve all bed items for a specified chromosome, from
    # this object. If there are no bed items for that chromosome (for this
    # track), then return an empty list:
    def get_BEDs(self, chrom):
        if (self.chroms.has_key(chrom)):
            return self.chroms[chrom]
        else:
            return []

    # Added on July 15th 2013:
    def getBEDsForRegion(self, region):
        """Returns a list of BED_line items whose start and end positions are
        located inside the specified region (specified as a BED_line object)."""

        chrom = region.getChrom()
        allBedItems = self.get_BEDs(chrom)

        # NOTE/FIXME: This assumes that the bed items in this track are of equal
        # size, and will return misleading output otherwise.
        # To fix this, I need to have separate indexes of the beds according
        # to chrom start and chrom end.

        # Get the index of the first bed item starting immediately inside the
        # specified region...

        idxOfBedClosestToStart = \
            list_funcs.find_closest_idx(region, allBedItems,
                                        get_disp_chrStartChrStart)
        bedClosestToStart = allBedItems[idxOfBedClosestToStart]
        if bedClosestToStart.getChromStart() > region.getChromStart():
            # That bed item is to the right of the start of the region.
            bedToGet_startIdx = idxOfBedClosestToStart
        else:
            # Will start from the next bed item over:
            bedToGet_startIdx = idxOfBedClosestToStart + 1

        # Get the index of the bed item ending just before the end of the
        # specified region...
        idxOfBedClosestToEnd = \
            list_funcs.find_closest_idx(region, allBedItems,
                                        get_disp_chrEndChrEnd )
        bedClosestToEnd = allBedItems[idxOfBedClosestToEnd]
        if bedClosestToEnd.getChromEnd() < region.getChromEnd():
            # That bed item is to the inside the region.
            bedToGet_endIdx = idxOfBedClosestToEnd
        else:
            # Will end at the previous bed item:
            bedToGet_endIdx = idxOfBedClosestToEnd - 1

        return allBedItems[bedToGet_startIdx:bedToGet_endIdx+1]

    def get_chroms(self):
        """Returns a list of the the chromosomes included in this bed object."""
        return self.chroms.keys()

    # Function to retrieve the item closest to the chromStart of the specified
    # BED line position, in log(N) time.
    def get_closest_item(self, query_item):
        # Retrieve the bed items for the specified chromosome:
        bed_items = self.get_BEDs(query_item.getChrom())

        # If no item is located on the same chromosome, then return None:
        if (len(bed_items) == 0):
            return None

        # Locate the item in the (sorted) list of bed items that is located
        # closest to the specified item:
        closest = list_funcs.find_closest(query_item, bed_items, get_disp)
        return closest

    # Function to retrieve the item whose chromEnd is closest to the
    # chromStart of the specified BED line position, in log(N) time.
    def get_closest_chromEnd(self, query_item):
        # Retrieve the bed items for the specified chromosome:
        bed_items = self.get_BEDs(query_item.getChrom())

        # If no item is located on the same chromosome, then return None:
        if (len(bed_items) == 0):
            return None

        # Locate the item in the (sorted) list of bed items that is located
        # closest to the specified item:
        closest = list_funcs.find_closest(query_item, bed_items, get_disp2)
        return closest

    # This function finds the closest bed item that occurs {upstream of and on
    # the same strand as} the specified query bed item:
    def get_closest_upstreamSameStrand(self, query_item):
        # Retrieve the bed items for the specified chromosome:
        bed_items = self.get_BEDs(query_item.getChrom())

        # If no item is located on the same chromosome, then return None:
        if (len(bed_items) == 0):
            return None

        # Iterate through the bed items from the *downstream* end of the
        # chromosome's beds, until the first {*upstream* bed on the same strand}
        # is encountered:
        if (query_item.getStrand() == "+"):
            # The query item is on the positive strand => Start scanning
            # from the *end* of the chromosome...
            bed_idx = len(bed_items) - 1

            # Logically, there must be at least one item in list, so the
            # following array indexing should not yield an error:
            curr_bed = bed_items[bed_idx]
            while ((bed_idx >= 0) and ((curr_bed.getStrand() == "-") or
                   (curr_bed.getChromStart() > query_item.getChromStart()))):
                curr_bed = bed_items[bed_idx]
                bed_idx = bed_idx - 1
        else:
            # The query item is on the negative strand => Start scanning
            # from the *end* of the chromosome...
            bed_idx = 0

            # Logically, there must be at least one item in list, so the
            # following array indexing should not yield an error:
            curr_bed = bed_items[bed_idx]
            while ((bed_idx < len(bed_items)) and
                   ((curr_bed.getStrand() == "+") or
                    (curr_bed.getChromEnd() < query_item.getChromEnd()))):
                curr_bed = bed_items[bed_idx]
                bed_idx = bed_idx + 1

        if ((bed_idx < 0) or (bed_idx >= len(bed_items))):
            return None
        else:
            return curr_bed


    def writeToFile(self, outFile):
        """Writes a bed file for this track."""
        for chrom in self.get_chroms():
            for bedItem in self.get_BEDs(chrom):
                outFile.write(bedItem.to_string() + "\n")
        outFile.flush()
                # FIXME: Currently, explicitly printing just the first 4 fields.
                # It would be nice to have a more flexible and general printing
                # mechanism.
#                outStr = \
#                    bedItem.getChrom() + " " + \
#                    bedItem.getChromStart() + " " + \
#                    bedItem.getChromEnd() + " " + \
#                    bedItem.getName() + "\n"
#                file.write(outStr)


    def toArray(self):
        """Returns a flattened version of this data structure; i.e. a list of
        BED_line objects."""
        arr = []
        for chrom in self.get_chroms():
            for bedItem in self.get_BEDs(chrom):
                arr.append(bedItem)
        return arr

    def sampleBEDs(self, nBEDs, sample="random"):
        # Sample the required number of bed items from all those in this
        # track...

        allBEDs = self.getAllBEDs()

        if (sample == "random"):
            # Sample items randomly:
            random.shuffle(allBEDs)
            return allBEDs[:nBEDs]
        elif (sample == "highScoreFilt"):
            # Select the highest scoring items after first ignoring the
            # top 5% of items:
            allBEDs.sort(compareScores, reverse=True)
            index5percent = int(len(allBEDs)/20.0)
            return allBEDs[:(index5percent+nBEDs)][-nBEDs:]
        else:
            # Select the highest scoring items:
            assert (sample == "highScore")
            allBEDs.sort(compareScores, reverse=True)
            return allBEDs[:nBEDs]


def gffPeaks2bed(gffFilename, outFilename):
    """Converts a specified gff file of ChIP-seq peaks into chip-seq peak-
    specific bed format."""

    # Introduced on October 1, 2010, for my chipseqQC analysis software:
    gffFile = open(gffFilename)
    outFile = open(outFilename, 'w')

    # FIXME: Should probably improve the parsing here in the future, to
    # introduce some error checking.

    # Generate an output line from each output line, and write it out:
    for line in gffFile:
        elems = line.split()
        # Minor error checking here. Just check that the line has the required
        # number of fields in it:
        if (len(elems) < 9):
            print >> sys.stderr, "Invalid input line in peaks gff file " + \
                gffFilename + ":", line
            sys.exit(1)

        chrom = "chr" + elems[0]
        chrStart = elems[1]
        chrEnd = elems[2]
        summit = elems[3]
        peakHeight = elems[4]
        peakVolume = elems[5]
        negCntrl = elems[6]
        fcRatio = elems[7]
        pVal = elems[8]

        extraInfo = peakVolume + "_" + \
            negCntrl + "_" + fcRatio + "_" + pVal
        # FIXME: Currently, I am excluding the extra info (above) for the peak,
        # rather than including it in the bed output line. Improve on this
        # later.
        outLine = chrom + " " + chrStart + " " + chrEnd + " NoName " + \
            peakHeight + " # " + pVal + " " + fcRatio + " " + summit + "\n"
        outFile.write(outLine)
    outFile.flush()
    outFile.close()


# This function returns the displacement between the chromStart of the first
# element, and the chromEnd of the second, assuming that they are located on
# the same chromosome as each other:
def get_disp2(bed1, bed2):
    assert(bed1.getChrom() == bed2.getChrom())
    return (bed1.getChromStart() - bed2.getChromEnd())


# This function returns the displacement from the chromStart of the second
# element to the chromStart of the first, assuming that they are located on
# the same chromosome as each other:
def get_disp_chrStartChrStart(bed1, bed2):
    assert(bed1.getChrom() == bed2.getChrom())
    return (bed1.getChromStart() - bed2.getChromStart())


# This function returns the displacement from the chromEnd of the second
# element to the chromEnd of the first, assuming that they are located on
# the same chromosome as each other:
def get_disp_chrEndChrEnd(bed1, bed2):
    assert(bed1.getChrom() == bed2.getChrom())
    return (bed1.getChromEnd() - bed2.getChromEnd())


# This function returns the displacement between the two elements, assuming
# that they are located on the same chromosome as each other:
def get_disp(bed1, bed2):
    if (bed1.getChromStart() > bed2.getChromEnd()):
        # Bed1 is "to the right" of bed2:
        return (bed1.getChromStart() - bed2.getChromEnd())
    elif (bed1.getChromEnd() < bed2.getChromStart()):
        # Bed1 is "to the left" of bed2 (so negative displacement):
        return -(bed2.getChromStart() - bed1.getChromEnd())
    else:
        # The two items must overlap => Report a displacement of zero:
        return 0


# This function chooses the bed line with the lowest score, out of the two
# options:
def get_lowest(line1, line2):
    if (line1.getScore() < line2.getScore()):
        return line1
    elif (line1.getScore() == line2.getScore()):
        print >> sys.stderr, "Lines have equal scores:\n", line1.to_string(), \
            "\n", line2.to_string()
        return line2
    else:
        return line2


# A comparator function. This function compares two bed lines objects according
# to their score values, and reports -1, 0, or 1, according to the outcome
# of the comparison.
def compareScores(line1, line2):
    if (line1.getScore() < line2.getScore()):
        return -1
    elif (line1.getScore() == line2.getScore()):
        return 0
    else:
        return 1


# A comparator function. Compares two chromosomal locations.
def chrom_compare(line1, line2):
    # Sort lexicographically according to chromosome first:
    if (line1.chrom > line2.chrom):
        return 1
    elif (line1.chrom < line2.chrom):
        return -1
    # If chromosomes are then same, then sort according to chromStart:
    elif (line1.chromStart > line2.chromStart):
        return 1
    elif (line1.chromStart == line2.chromStart):
        return 0
    else:
        return -1


class BED_line(object):
    """A line of data in BED format. Modified on 6th of October 2010: Introduced
    a comments field. Now, any non-compulsory field will be regarded as a
    comment if it occurs after a "#" character. I did this to allow me
    to extend the class and to include extra information in the comments field
    for the extending class peak_BED_line."""

    def __init__(self, line, track=None):
        items = line.split()
        # Introduced on 29th October 2010 to allow storage of wiggle track data
        # for an individual BED item:

        # August 1st 2013:
        # Removing reference to parent track; I suspect it could create a
        # "reference cycle", which would make garbage collection not happen:
        #self.track=track # The track to which this BED line belongs (if any)

        self.wigVals = {}
        try:
            # Parse the compulsory fields (chrom, start, end):
            self.chrom = items[0]
            self.chromStart = int(items[1])
            self.chromEnd = int(items[2])

            # Parse the remaining fields on the line...
            self.comment = None
            self.name = None
            self.score = None
            self.strand = None

            # For each remaining item from the 4th onwards, check if it
            # starts with "#". If so, store the remaining sections of the string
            # as a comment and stop parsing the line. Otherwise, store the
            # field as the appropriate attribute in the class:
            item_idx = 3
            while ((self.comment == None) and (len(items) > item_idx)):
                currItem = items[item_idx]
                if (currItem[0] == "#"):
                    # The remaining fields on the line comprise a comment.
                    # Reconstitute the string from those fields, and store
                    # the resulting comment:
                    commentStr = reduce(lambda tok1, tok2: tok1 + " " + tok2,
                                        items[item_idx:])[1:]
                    self.comment = commentStr
                else:
                    # The field is a regular optional bed line field. Parse
                    # it as such:
                    if (item_idx == 3):
                        self.name = currItem
                    if (item_idx == 4):
                        self.score = float(currItem)
                    if (item_idx == 5):
                        if (currItem == "None"):
                            # Added on 15th October 2010: Allow "None" as
                            # a valid strand string:
                            self.strand = None
                        else:
                            if (currItem != '-' and currItem != '+'):
                                raise ValueError("Invalid strand string: " +
                                                 currItem)
                        self.strand = currItem
                    
                item_idx = item_idx + 1

        except ValueError:
            # If the BED line is not formatted correctly, then raise the
            # exception to the calling code.
            raise

            # FIXME: The following fields of BED objects are not being
            # used in my current applications. => Comment them out (ie don't
            # make them attributes of the class) for the time being.
#             if (len(items) > 6):
#                 self.thickStart = items[6]
#             else:
#                 self.thickStart = None
#             if (len(items) > 7):
#                 self.thickEnd = items[7]
#             else:
#                 self.thickEnd = None
#             if (len(items) > 8):
#                 self.itemRgb = items[8]
#             else:
#                 self.itemRgb = None
#             if (len(items) > 9):
#                 self.blockCount = items[9]
#             else:
#                 self.blockCount = None
#             if (len(items) > 10):
#                 self.blockSizes = items[10]
#             else:
#                 self.blockSizes = None
#             if (len(items) > 11):
#                 self.blockStarts = items[11]
#             else:
#                 self.blockStarts = None

    # Converts this bed line to a string:
    def to_string(self):
        str_rep = self.chrom + " " + str(self.chromStart) + " " + \
            str(self.chromEnd) + " " + str(self.name) + " " + str(self.score) \
            + " " + str(self.strand)
        return  str_rep

    # August 1st 2013:
    # Removing reference to parent track; I suspect it could create a
    # "reference cycle", which would make garbage collection not happen:
    #def getTrack(self):
    #    return self.track

    def getChrom(self):
        return self.chrom

    def getChromStart(self):
        return self.chromStart

    def setChromStart(self, newStart):
        self.chromStart = newStart

    def setChromEnd(self, newEnd):
        self.chromEnd = newEnd

    def getChromEnd(self):
        return self.chromEnd

    def getLen(self):
        """Return the length of this genomic interval."""
        return self.getChromEnd() - self.getChromStart() + 1

    def getScore(self):
        return self.score

    def setScore(self, score):
        self.score = score

    def setName(self, name):
        self.name = name

    def getStrand(self):
        return self.strand

    def getName(self):
        return self.name

    # Returns true if this bed line overlaps the specified item's chromStart.
    # Returns false otherwise:
    def overlaps(self, bed_item):
        if (self.getChrom() != bed_item.getChrom()):
            return False
        start = self.getChromStart()
        end = self.getChromEnd()
        if ((start <= bed_item.getChromStart()) and \
            (end >= bed_item.getChromStart())):
            return True
        else:
            return False

    def overlaps2(self, bed_item):
        """Returns true if this bed line overlaps any of the specified item's
        region. Returns false otherwise:"""
        if (self.getChrom() != bed_item.getChrom()):
            return False
        start = self.getChromStart()
        end = self.getChromEnd()
        if ((start <= bed_item.getChromStart()) and
            (end >= bed_item.getChromStart())):
            return True
        elif ((start <= bed_item.getChromEnd()) and
              (end >= bed_item.getChromEnd())):
            return True
        elif ((start >= bed_item.getChromStart()) and
              (end <= bed_item.getChromEnd())):
            return True
        else:
            return False

    def setWig(self, wigKey, wigArray):
        """Sets wiggle data for this BED object. Stores the specified array
        under the specified object key, which will be used as a dictionary
        key."""
        self.wigVals[wigKey] = wigArray

    def getWig(self, wigKey):
        """Returns the wiggle data for the specified key."""
        return self.wigVals[wigKey]

    def spansBED(self, bedTrack):
        """"Returns true if there is at least one {bed item in the specified
        bedTrack} whose end occurs within the bounds of this bed item."""
        regStart = self.getChromStart()
        regEnd = self.getChromEnd()
        regMiddle = int((regStart + regEnd)/2.0)
        tmpBedLineStr = self.getChrom() + " " + str(regMiddle) + " " + \
            str(regMiddle)
        halfRegSize = (regEnd - regStart)/2.0
        regMiddle_BED = BED_line(tmpBedLineStr)
        bedMiddleClosestPeak = bedTrack.get_closest_item(regMiddle_BED)
        if ((bedMiddleClosestPeak != None) and \
            (abs(get_disp(bedMiddleClosestPeak, regMiddle_BED)) < \
             halfRegSize)):
            return True
        else:
            return False

    def spanningBED(self, bedTrack):
        """"Returns the closest item in bedTrack if there is at least one {bed
        item in the specified bedTrack} that spans the middle of this bed
        item. Otherwise returns None"""
        regStart = self.getChromStart()
        regEnd = self.getChromEnd()
        regMiddle = int((regStart + regEnd)/2.0)
        tmpBedLineStr = self.getChrom() + " " + str(regMiddle) + " " + \
            str(regMiddle)
        regMiddle_BED = BED_line(tmpBedLineStr)
        bedMiddleClosestItem = bedTrack.get_closest_item(regMiddle_BED)
        if ((bedMiddleClosestItem != None) and
            (bedMiddleClosestItem.getChromStart() < regMiddle) and
            (bedMiddleClosestItem.getChromEnd() > regMiddle)):
            return bedMiddleClosestItem
        else:
            return None

    def getBedCoordStr(self):
        """Returns a single-word representation of the coordinates of this
        bed object."""

        return str(self.getChrom()) + ":" + str(self.getChromStart()) + "-" + \
            str(self.getChromEnd())

    # NOTE: Added this functionality on 1st October 2010 in order to support
    # the chipseq QC analysis code:
    # DEPRECATED ON 14TH OCTOBER 2010, as I have refactored my code =>
    # Commenting this out.
#     def seq2fasta(self, genomeFilename, outFilename, noZero=False):
#         """This method retrieves the genomic sequence corresponding to the
#         BED line's coordinates, and writes that sequence out to the specified
#         file, in fasta format."""

#         # noZero is a boolean specifying whether a zero-length sequence is
#         # to be output.

#         # FIXME: Currently, I am assuming the input reference genome chromosome
#         # sequences are named just by number, without a preceding "chr". This
#         # assumption is very volatile, and I should improve the mechanism
#         # here to make it robust in the future - e.g. by parsing part of
#         # the genome fasta file to ensure it follows that format.
#         chrLocStr = self.getChrom()[3:] + ":" + str(self.getChromStart()) + \
#             "-" + str(self.getChromEnd())

# #        print >> sys.stderr, \
# #            "seq2fasta called. Attempting to run pysam.faidx for the inputs", \
# #            genomeFilename, "and", chrLocStr

#         seqGetter = getSeqFa.seqRetriever(genomeFilename)
#         seq = seqGetter.getSequence(chrLocStr)

#         # Output the sequence in fasta format...
#         seqNameLine = ">" + chrLocStr + "\n"
#         outFile = open(outFilename, 'a')

#         if (noZero):
#             if (len(seq) > 0):
#                 # Sequence is non-zero length => can print it:
#                 outFile.write(seqNameLine)
#                 outFile.write(seq + "\n")
#         else:
#             # Output the sequence, irrespective of it's length:
#             outFile.write(seq + "\n")
#             outFile.write(line)

#         outFile.flush()
#         outFile.close()


def peaksFcCmp(peakBED1, peakBED2):
    """Comparator function that compares two peak beds on fold-change."""
    if (peakBED1.foldChange < peakBED2.foldChange):
        return -1
    elif (peakBED1.foldChange == peakBED2.foldChange):
        return 0
    else:
        return 1


class peak_BED_line(BED_line):
    """A class representing a special kind of BED line, which represents a
    chip-seq peak. This line has three extra attributes (pValue, foldChange,
    and summit), which are stored as space-delimited items in that order in
    the comment string of the input line. I have done it this way because
    the standard bed line does not include these attributes as defined by
    UCSC."""

    def __init__(self, line, track=None):
        # Instantiate all the fields on the line:
        super(peak_BED_line, self).__init__(line, track=track)

        # Parse the comment string of the parent class, and extract the 3
        # required fields from it. Report an error if the comment string
        # formatting is invalid...
        try:
            # Check that the comments field has been instantiated and that
            # it has three items. If it doesn't, report this input error as
            # an exception:
            if ((self.comment == None) or (len(self.comment.split()) != 3)):
                raise ValueError(\
                    "Invalid input line for peak_BED_line creation:" + \
                        str(self.comment))

            toks = self.comment.split()
            self.pVal = float(toks[0])
            self.foldChange = float(toks[1])
            self.summit = int(toks[2])

        # Parse the comment line, storing 
        except ValueError:
            # If the BED line is not formatted correctly, then raise the
            # exception to the calling code.
            raise

        # The peak sequence is also introduced as an attribute of this
        # BED_line subclass:
        peakSequence = None

    def get_pVal(self):
        return self.pVal

    def get_fc(self):
        return self.foldChange

    def get_summit(self):
        return self.summit

    def to_string(self):
        """Returns a string representation of this peak bed object. Simply
        appends the additional variables onto the end of the string generated
        by the superclass."""
        superStrRep = super(peak_BED_line, self).to_string()
        str_rep = superStrRep + " # " + str(self.get_pVal()) + " " + \
            str(self.get_fc()) + " " + str(self.get_summit())
        return  str_rep

    # August 1st 2013:
    # Removing reference to parent track; I suspect it could create a
    # "reference cycle", which would make garbage collection not happen:
    # def getRefGenomeFilename(self):
    #     return self.getTrack().getRefGenomeFilename()

    # September 18th 2013:
    # Re-adding the getPeakSequence method, but this time the user has
    # to pass in a seqGetter object:

    def getPeakSequence(self, flankWidth, seqGetter):
        """Returns the peak's sequence, corresponding to the sequence
        surrounding the summit with the specified flanking sequence width."""

        # Generate a new chromosomal loc string, centred on the peak's
        # summit, and with the specified width:
        startPos = str(self.get_summit() - flankWidth)
        endPos = str(self.get_summit() + flankWidth)
        if (int(startPos) < 0):
            print >> sys.stderr, "WARNING: Ignoring peak as extended start " + \
                "is below 0:", self.to_string()
            return None

        centredPeakLoc = self.getChrom() + ":" + startPos + "-" + endPos

        # Retrieve the genomic sequence corresponding to that location,
        # and set this peak's sequence accordingly...

        # March 4th 2013: Trim the start of the chromosome string name
        # if needed. Figure out whether to trim based on the genome
        # filename contents. I don't like this solution but it is the
        # best solution I can come up with in the current design, since
        # option passing is not flexible enough in chipseq_analysis.py:
        firstGenomeLine = open(seqGetter.seqFilename).readline()
        chrLocQueryStr = None
        if firstGenomeLine[:4] == ">chr":
            chrLocQueryStr = centredPeakLoc
        else:
            chrLocQueryStr = centredPeakLoc[3:]
        try:
            peakSequence = seqGetter.getSequence(chrLocQueryStr)
            # Ensure that the peak sequence is the correct size given
            # the flanking width specified. If it is not, then set the
            # peakSequence to null and report this as a warning:

            requiredSeqSize = (flankWidth*2)+1
            if (len(peakSequence) != requiredSeqSize):
                print >> sys.stderr, "WARNING!: Genomic sequence for " + \
                    "peak " + self.to_string() + " is not the required size" + \
                    " of " + str(requiredSeqSize) + ". Setting sequence to null" + \
                    " instead."
                return None
            else:
                return peakSequence

        except ValueError, e:
            # If the peak's sequence could not be retrieved, report this
            # error (a possible user input error) but do not crash:
            print >> sys.stderr, "WARNING!: Could not retrieve sequence" + \
                " for peak " + self.to_string() + "."
            
            # Sequence retrieval failed for this peak for some reason =>
            # Make sure that the peakSequence is returned as null
            # to flag this fact:
            return None

    def writeFasta(self, outFile, flankWidth, seqGetter):
        """Write the sequence for this peak to the specified output fasta
        file. The sequence name written is exactly the same as the string
        returned by the parent class's method getBedCoordStr()."""
        seqName = ">" + self.getBedCoordStr()

        # NOTE: Potentially confusing point: The sequence returned probably
        # does *not* correspond to the printed sequence name; This is because
        # "getPeakSequence" retrieves the region flanking the peak summit,
        # whereas the sequence name specifies the start and end of the peak
        # itself.

        seq = self.getPeakSequence(flankWidth, seqGetter)
        outFile.write(seqName + "\n")
        outFile.write(seq + "\n")

    def copy(self):
        """Introduced on January 21st 2011 to facilitate spamo analysis.
        Returns a copy of this peak_BED_line object."""

        # Generate a string representing the line, which is required by
        # the constructor of this class:
        strRep = self.to_string()

        # Generate the copy and return it:
        copy = peak_BED_line(strRep)
        copy.track = self.getTrack()
        return copy

    def setStart(self, newStart):
        """Introduced on January 21st 2011 to facilitate spamo analysis.
        Sets the start coordinates of this BED_line object. NOTE: This
        is a bit dodgy as I had originally planned for the start and end
        of a given BED_line object to be immutable."""
        self.chromStart = newStart

    def setEnd(self, newEnd):
        """Introduced on January 21st 2011 to facilitate spamo analysis.
        Sets the end coordinates of this BED_line object. NOTE: This
        is a bit dodgy as I had originally planned for the start and end
        of a given BED_line object to be immutable."""
        self.chromEnd = newEnd


class FixedStep_Track_From_BED:
    """A fixedStep wiggle track class that is created by reading in a bed-
    formatted file."""
    def __init__(self, file, bg_score):
        # A hashtable with chromosomal locs as keys and scores from the
        # bed items as values. Supports O(1) retrieval of values on the basis
        # of chromosomal location, assuming that no collisions occur.
        self.values = {}

        # The background score:
        self.bg_score = bg_score

        # Parse the contents of the specified file to produce the wiggle
        # track object:
        try:
            for line in file:
                elems = line.split()

                # Will just ignore the line if it is empty; empty lines are
                # valid:
                if (len(elems) > 0):
                    # Only add values for the line if it is a chromosomal
                    # location:
                    first_elem = elems[0]
                    if (first_elem[0:3] == "chr"):
                        # Each bed entry must have a score entry if the bed
                        # object is going to be converted into a fixedStep
                        # wiggle track object:
                        if (len(elems) < 5):
                            raise ValueError(
                                "Line invalid for FixedStep_Track_From_BED:\n" \
                                + line)
                        chrom = elems[0]
                        start = int(elems[1])
                        end = int(elems[2])
                        score = float(elems[4])
                        # Generate a list of chromosomal locations, spanning the
                        # specified start to the specified end:
                        for loc in (range(start,(end - 1))):
                            curr_pos = (chrom, loc)
                            # Add each resulting chromosomal location string to
                            # the values hashtable, with the items score as the
                            # value:
                            self.values[curr_pos] = score
                
        except IOError:
            # If the file does not open, throw this exception upwards, rather
            # than dealing with it here:
            raise

        except ValueError:
            # A value error here means that a BED_line was invalid.
            raise


    # Retrieves the wiggle track's value at the specified chromosomal position.
    def get_value(self, chrom_name, pos):
        if (self.values.has_key((chrom_name,pos))):
            return self.values[(chrom_name,pos)]
        else:
            return self.bg_score # No explicit value => Return background value.


class FixedStep_Track:
    """A fixedStep wiggle track."""

    # Added max_chrom_size on 24th July 2008. It specifies the maximum size
    # a chromosome can be. If it is specified, then the array of float for
    # each chromosome is pre-allocated and filled in. If left as null, the
    # array is grown using append, instead.
    def __init__(self, file, max_chrom_size=None):
        print >> sys.stderr, "Reading file: ", file.name # Trace
#         # Insist that the first line must define the fixedStep interval. Read
#         # in that info:
#         first_line = file.readline()
#         elems = first_line.split()
#         if (elems[0] != "Track"):
#             raise ValueError("First line of file " + file.name + \
#                              " must start with \"track\"." + items[5])
        # The arrays of data, corresponding to the different chromosomes, will
        # be stored in a dictionary for easy retrieval:
        self.chroms = {}

        # Discard the lines at the start of the file until a "fixedStep" line
        # is found:
        curr_line = file.readline()
        while ((curr_line != "") and (len(curr_line.split()) > 0) and \
               (curr_line.split()[0] != "fixedStep")):
            curr_line = file.readline()

        # Keep reading data and placing it in "Chromosome" objects until no
        # more lines are available in the file:
        while (curr_line != ""):
            # Current line is a "fixedStep" definition line. Retrieve the
            # chrom, start, and step values from it:
            elems = curr_line.split()
            chrom_str = elems[1]
            start_str = elems[2]
            step_str = elems[3]
            chrom = chrom_str.split("=")[1]
            start = int(start_str.split("=")[1])
            step = int(step_str.split("=")[1])

            # Read in data until the end of the file is reached or a new
            # chromosome is reached:
            if (max_chrom_size == None):
                chrom_data = array.array('f', [])
                print >> sys.stderr, "STARTING WITH ARRAY OF SIZE ZERO."
            else:
                # 25th July 2008: Forced to use python arrays, to reduce
                # memory usage.
                
                # The maximum size for the chromosome has been specified:
                # Default value; -100 value seems OK...
                default_val = -100
                chrom_data = array.array('f', [default_val]*max_chrom_size)
                print >> sys.stderr, "STARTING WITH EMPTY ARRAY OF SIZE", \
                    max_chrom_size
                curr_chrom_idx = 0
            curr_line = file.readline()
            while ((curr_line != "") and (len(curr_line.split()) > 0) and \
                   (curr_line.split()[0] != "fixedStep")):
                if (max_chrom_size == None):
                    chrom_data.append(float(curr_line))
                else:
                    chrom_data[curr_chrom_idx] = float(curr_line)
                    curr_chrom_idx = curr_chrom_idx + 1
                curr_line = file.readline()
            print >> sys.stderr, "FINISHED PROCESSING FIXEDSTEP BLOCK."

            # Create a new chromosome, containing the values. These wiggle
            # tracks are only allowed to have one segment of data per
            # chromosome, for simplicity. Note that hence this class
            # does not fully implement a wiggle track. Check this here:
            if (self.chroms.has_key(chrom)):
                raise ValueError( \
                    "More than one segment of data for a single chromosome: " \
                    + chrom + " . File: " + file.name)
            self.chroms[chrom] = Chromosome(chrom, start, step, chrom_data)


    # Retrieves the wiggle track's value at the specified chromosomal position.
    # Returns None if the wiggle track does not record a value at that location:
    def get_value(self, chrom_name, pos):
        # Retrieve the specified chromosome:
        curr_chrom = self.chroms[chrom_name]

        return curr_chrom.get_value(pos)


    # Prints the data for the specified chromosome, for this fixedStep Wiggle
    # track:
    def print_chrom(self, outfile, chrom_name, chrom_len):
        assert(self.chroms.has_key(chrom_name))
        # The length of the chromosome must be specified, as this is not
        # specified when the track is first created. FIXME: Could improve that.
        print >> outfile, self.chroms[chrom_name].to_string(chrom_len)
        

class Chromosome:
    def __init__(self, chrom, start, step, chrom_data):
        self.chrom = chrom
        self.start = start # Start of the data values in this chromosome
        self.step = step
        self.data = chrom_data


    def get_value(self, pos):
        # Calculate index of specified position into the array:
        arr_idx = int(math.floor(float(pos - self.start)/self.step))
        # int(math.ceil(float(pos - self.start)/self.step) - 1)

        # If there is an item of data at that location, return it. Otherwise
        # return None to indicate that the location did not have a value
        # in this wiggle track:
        if ((arr_idx >= 0) and (arr_idx < len(self.data))):
            return self.data[arr_idx]
        else:
            return None


    def to_string(self, chrom_len):
        out_str = "CHROM: " + str(self.chrom) + "\n"
        out_str = out_str + "START: " + str(self.start) + "\n"
        out_str = out_str + "STEP: " + str(self.step) + "\n"
        for pos in range(1, chrom_len + 1):
            out_str = out_str + str(self.get_value(pos)) + "\n"
        return out_str


def file_is_fixedStep(filename):
    """Scans the file, and returns if it is in fixedStep format, false
    otherwise"""
    
    file = open(filename)

    for line in file:
        elems = line.split()
        if (elems[0] == "fixedStep"):
            return True

    return False


def file_is_bed(filename):
    """Scans the file, and returns if it is in BED format, false
    otherwise"""
    
    file = open(filename)

    for line in file:
        elems = line.split()
        if (elems[0][0:3] == "chr"):
            return True

    return False


def getChromSizes(infoFile):
    """Returns a hashtable containing the lengths of the chromosomes, as
    recorded in the specified file."""
    sizeHash = {}
    for line in infoFile:
        elems = line.split()
        sizeHash[elems[0]] = int(elems[1])

    return sizeHash
