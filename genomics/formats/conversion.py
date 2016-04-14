#!/usr/bin/env python

import genomics.formats.bed as bed
import cmdlineProgs, pdb

"""A module for performing various genomics-oriented file format
conversions. Interfaces with packages such as samtools and IGV tools."""


def macs2bed(inputXlsFilename, outputBedFilename):
    """Converts the specified MACS-formatted XLS file into BED format.
    Overwrites the specified bed output file."""

    bedOutfile = open(outputBedFilename, 'w')
    for line in open(inputXlsFilename).xreadlines():
        peakIdx = 1
        if ((len(line) > 3) and (line[0] != "#") and (line[:3] != "chr")):
            elems = line.strip().split("\t")
            # A bit ugly: Combining all the summit, q-value, peak size etc.
            # info into the name string:
            peakName = "peak_" + str(peakIdx) + "#" + \
                '_'.join([elems[3], elems[4], elems[6], elems[7], elems[8], elems[-1]])
            bedLine = " ".join([elems[0], elems[1], elems[2], peakName, elems[5]])
            print >> bedOutfile, bedLine
            peakIdx += 1

    bedOutfile.close()


def macs2bedTrack(inputXlsFile):
    """Parses the specified MACS xls file and returns a bed.BED_Track object
    with peakBED items."""

    peaksList = []
    for line in inputXlsFile.xreadlines():
        peakIdx = 1
        if ((len(line) > 3) and (line[0] != "#")):
            elems = line.strip().split("\t")

            if elems[0] != "chr":
                # Generate a "peak_BED_line" object representing this peak:
                summit = str(int(elems[1]) + int(elems[4]))
                currPeakStr = " ".join([elems[0], elems[1], elems[2], "#", elems[-1], elems[-2], summit])
                
                currPeak = bed.peak_BED_line(currPeakStr)
                peaksList.append(currPeak)

    peaksTrack = bed.BED_Track(peaksList, peaksTrack=True)
    return peaksTrack


def wig2tdf(wigFileloc, tdfFileloc, genome = "hg19", verbose="False"):
    """Converts the specified wig file to tdf file."""

    # Run igvtools "toTDF" and "index" to accomplish this...
    cmdStr = "/mnt/oldHome/tom/Work/externalProgs/IGVTools/igvtools toTDF -f min,max,mean " + wigFileloc + " " + tdfFileloc + " " + genome
    cmdlineProgs.runCommand(cmdStr, verbose)


def bam2bed(bamFilename, bedFilename, regionStr="", verbose=False):
    """This currently uses samtools and bedtools. Region is a string
    specifying chromosomal coordinates of the region to retrieve reads from."""

    cmdStr = "samtools view -b " + bamFilename + " " + regionStr + \
        " | bedtools bamtobed -i  - > " + bedFilename
    cmdlineProgs.runCommand(cmdStr, verbose)


def bed2bam(bedFilename, bamFilename, genomeFile, verbose=False):
    """This currently uses bedtools."""

    cmdStr = "bedtools bedtobam -i " + bedFilename + " -g " + genomeFile + \
        " > " + bamFilename
    cmdlineProgs.runCommand(cmdStr, verbose)
