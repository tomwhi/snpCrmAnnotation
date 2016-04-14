#!/usr/bin/env python

import commands, MySQLdb, gzip, pdb, sys
import utility, cmdlineProgs


def extractVCFsubsetTupleInput(outputVcfGzFilename, querySNPsTups, allSNPsVCFfileloc):
    """The same as extractVCFsubset(), except that it accepts a list
    of SNP tuples as inputs. Expands their locations by 1bp to the left but
    still only returns the exact SNPs requested."""

    # Generate the temporary BED file:
    tmpBEDfilename = utility.makeTempFilename("tmpBED")
    tmpBEDfile = open(tmpBEDfilename, 'w')
    for tup in querySNPsTups:
        chrom = tup[0]
        pos = int(tup[1])
        name = tup[2]
        print >> tmpBEDfile, chrom, pos-1, pos, name
    tmpBEDfile.close()

    # Extract the VCF entries:
    extractVCFsubset(outputVcfGzFilename, tmpBEDfilename, allSNPsVCFfileloc)

    cmdlineProgs.deleteFiles([tmpBEDfilename])


def extractVCFsubset(outputVcfGzFilename, querySNPsbedFileloc, allSNPsVCFfileloc):
    """Extracts the specified dictionary of SNPs ((keys are (chrom/pos/name)
    tuples) from the specified .vcf.gz tabixed file of all SNPs, and outputs to
    the tabixed data to the specified vcf.gz file.

    The input bed file and input vcf file must have the same chromosome name
    naming convention (i.e. whether to include "chr" at the start or not)."""

    # Use the header from the file of all SNPs:
    tmpVCFfilename = utility.makeTempFilename("tmpVCF")
    cmdStr = "zcat " + allSNPsVCFfileloc + " | head -10000 | awk '$1 ~ /#/' > " + tmpVCFfilename
    cmdlineProgs.runCommand(cmdStr, True)
    tmpVCFfile = open(tmpVCFfilename, 'a')

    # Use tabix to retrieve the relevant data. Need to do this
    # in two steps: first retrieve with tabix, and then filter the resulting
    # text VCF file:
    tmpVCFlinesFilename = utility.makeTempFilename("tmpVCFlines")
    cmdStr2 = "tabix -B " + allSNPsVCFfileloc + " " + querySNPsbedFileloc + " > " + tmpVCFlinesFilename
    cmdlineProgs.runCommand(cmdStr2, True)

    querySNPdict = {}
    for line in open(querySNPsbedFileloc).xreadlines():
        elems = line.strip().split()
        currTup = (elems[0], elems[2], elems[3])
        querySNPdict[currTup] = 1

    for vcfLine in open(tmpVCFlinesFilename).xreadlines():
        elems = vcfLine.strip().split()
        lineAsTup = (elems[0], elems[1], elems[2])
        if querySNPdict.has_key(lineAsTup):
            print >> tmpVCFfile, vcfLine.strip()

    tmpVCFfile.close()

    cmdStr3 = "cat " + tmpVCFfilename + " | bgzip > " + outputVcfGzFilename
    cmdlineProgs.runCommand(cmdStr3, True)

    cmdStr4 = "tabix -p vcf " + outputVcfGzFilename
    cmdlineProgs.runCommand(cmdStr4, True)
    
    cmdlineProgs.deleteFiles([tmpVCFfilename])

    # Hacky: Cope with differences in chromosome names ("chr" or not) here:
    #firstSNPchrom = snpDict.keys()[0][0]
    #queryStartsWithChr = True
    #if (len(firstSNPchrom) < 4) or (firstSNPchrom[:3] != "chr"):
    #    queryStartsWithChr = False

    # Problem: Need to initiate a new pipe from the specified all-snps gzipped
    # file to a python process to read all the snps. For the sake of cleanliness,
    # I'd prefer that to be this python process. Current solution: I'll use
    # python's zcat functionalities:
    # allSNPsFile = gzip.open(allSNPsVCFfileloc)
    # pdb.set_trace()
    # currLine = allSNPsFile.readline()
    # while currLine != "":
    #     elems = currLine.strip().split()
    #     targetChrom = elems[0]
    #     targetChromForComparison = targetChrom
    #     targetStartsWithChr = True
    #     if (len(targetchrom) < 4) or (targetchrom[:3] != "chr"):
    #         targetStartsWithChr = False
    #     if (queryStartsWithChr and not targetStartsWithChr):
    #         targetChromForComparison = "chr" + targetChrom
    #     elif (not queryStartsWithChr and targetStartsWithChr):
    #         targetChromForComparison = targetChrom[3:]
    #     pos = elems[1]
    #     name = elems[2]
    #     if snpDict.has_key((targetChromForComparison, pos, name)):
    #         print >> tmpVCFfile, currLine.strip()
    #     currLine = allSNPsFile.readline()

    # tmpVCFfile.close()


def vcfFile2chromPosNameIndex(vcfFile):
    """Produces a dictionary linking (chrom, pos, name) tuples to'
    vcffile lines."""

    chromPosName2vcfLine = {}
    for line in vcfFile.readlines():
        if len(line.strip()) > 0 and line[0] != "#":
            elems = line.strip().split()
            chromPosName2vcfLine[(elems[0], int(elems[1]), elems[2])] = " ".join(elems[:5])

    return chromPosName2vcfLine


def getSNV(chrom, refPos, name, chromPosName2vcfLine):
    """Retrieves a specified SNV from a specified vcf file and returns
    it as an SNV object."""

    # Currently just iteracte therough the vcf file until the desired
    # snv is found...
    if not chromPosName2vcfLine.has_key((chrom,refPos,name)):
        pdb.set_trace()
        x = 1
    assert chromPosName2vcfLine.has_key((chrom,refPos,name))
    vcfLine = chromPosName2vcfLine[(chrom,refPos,name)]
    elems = vcfLine.strip().split()
    alleleSeqs = [elems[3]] + elems[4].split(",")

    snv = SNV(chrom, int(refPos), name, alleleSeqs)
    return snv


class SNV:
    def __init__(self, chrom, refPos, name, alleleSeqs):
        self.chrom = chrom
        self.refPos = refPos
        self.name = name

        # The input array must be >= 2 long. The first element is
        # the reference genome sequence:
        assert(len(alleleSeqs) >= 2)
        self.alleleSeqs = alleleSeqs

    def getChrom(self):
        return self.chrom

    def getRefPos(self):
        return self.refPos

    def getAlleleSeq(self, alleleIdx):
        """Returns the allele sequence for the specified allele,
        with respect to the genomic plus strand."""
        return self.alleleSeqs[alleleIdx]

    def getRefAlleleSeq(self):
        return self.getAlleleSeq(0)

    def getAlleleID(self, alleleIdx):
        return self.getLocStr() + "_" + self.name + "_" + self.getAlleleSeq(alleleIdx)

    def getAlleleIDs(self):
        return map(lambda alleleIdx: self.getAlleleID(alleleIdx), self.getAlleleIdxs())

    def getAlleleIdxs(self):
        return range(len(self.alleleSeqs))

    def getLocTup(self):
        return (self.chrom, self.refPos)

    def getSNPtup(self):
        return (self.chrom, self.refPos, self.name)

    def getLocStr(self):
        return self.chrom + "_" + str(self.refPos)

    def toString(self):
        return self.chrom + "_" + str(self.refPos) + "_" + self.name

    def toBedString(self):
        """Generate a bed item representing this SNV's location."""

        # HACK: Have to use location as name, in order to be able to parse
        # locations from the output of bigWigAverageOverBed. How annoying:
        locName = self.chrom + "_" + str(self.refPos) + "_" + str(self.refPos + 1)

        return self.chrom + "\t" + str(self.refPos) + "\t" + str(self.refPos + 1) + "\t" + locName

    def getName(self):
        return self.name

    def getAlleleWindowSeq(self, alleleIdx, seqGetter, width):
        """Returns the plus strand genomic sequence that the specified
        allele would have, with a width-1 overhang on either side."""

        # FIXME: Could probably do this more efficiently and save
        # half the time. Not currently in an inner loop so not bothering.

        assert(width >= 1)

        # Get the width of the reference allele:
        refAlleleWidth = len(self.getAlleleSeq(0))

        # Get the genomic sequences flanking the reference allele:
        flankSize = width - 1
        leftFlankPos = self.chrom + ":" + str(self.refPos - flankSize) + "-" + \
            str(self.refPos - 1)
        leftFlankStr = seqGetter.getSequence(leftFlankPos)
        # June 20th 2013: Fixing off-by-one bug:
        rightFlankPos = self.chrom + ":" + str(self.refPos + refAlleleWidth) + \
            "-" + str(self.refPos + refAlleleWidth + flankSize - 1)
        rightFlankStr = seqGetter.getSequence(rightFlankPos)
        
        # October 2 2014: Make sure sequence is lower case except for SNP nucleotides,
        # which are upper case:
        return leftFlankStr.lower() + self.getAlleleSeq(alleleIdx).upper() + rightFlankStr.lower()

    def writeAlleleWindowSeqs(self, outFile, seqGetter, width):
        """Writes out the allele sequences for this SNV to the specified
        output stream, with width-1 overhang on either side. seqGetter
        must be for the same reference genome as this SNV is wrt."""

        # Write out each of the alternative alleles:
        for alleleIdx in self.getAlleleIdxs():
            currAlleleID = self.getAlleleID(alleleIdx)
            print >> outFile, ">" + currAlleleID
            print >> outFile, self.getAlleleWindowSeq(alleleIdx, seqGetter, width)

    def getLDinfo(self, options):
        """Gets an array of tuples, of the form (ldStatistic, SNV), for each SNV
        that is in LD with this SNV. Retrieves this information from the
        snpEffect mysql database whose login info is specified in the options
        input."""

        r2thresh = float(options.r2thresh)
        snpEffectDb=MySQLdb.connect(db=options.snpEffectDb,
                                    user=options.uname,
                                    read_default_file="~/.my.cnf")

        # Get all SNVs associated with this SNV by
        # LD, using queries to the mysql database...
        ldSNPtups = []

        cursor=snpEffectDb.cursor()
        cmdStr = "drop table if exists currLDsnps;"
        cursor.execute(cmdStr)

        cmdStr = "create table currLDsnps select chrom, snp2pos as pos, r2 from snpLD where snp1pos = \"" + str(self.refPos) + "\" and chrom = \"" + self.chrom + "\";"
        cursor.execute(cmdStr)

        cmdStr = "insert into currLDsnps select chrom, snp1pos as pos, r2 from snpLD where snp2pos = \"" + str(self.refPos) + "\" and chrom = \"" + self.chrom + "\";"
        cursor.execute(cmdStr)

        cmdStr = "SELECT 1000genomes.*, currLDsnps.r2 from 1000genomes inner join currLDsnps on 1000genomes.chrom = currLDsnps.chrom and 1000genomes.pos = currLDsnps.pos;"
        cursor.execute(cmdStr)
        rows = cursor.fetchall()

        ldSNPtups = ldSNPtups + list(rows)

        # Filter on r2 threshold and convert tuples to SNV objects:
        ldSNPtups_filt = []
        for snvTup in ldSNPtups:
            alleleSeqs = [snvTup[3]] + snvTup[4].split(",")
            currSNV = SNV(snvTup[0], int(snvTup[1]), snvTup[2], alleleSeqs)
            snvTup = (currSNV, snvTup[-1])
            if snvTup[-1] > r2thresh:
                ldSNPtups_filt.append(snvTup)
        
        ldSNPtups_filt.append((self, -1))

        return ldSNPtups_filt

    def toVCFline(self):
        altAlleleSeqsStr = reduce(lambda tok1, tok2: tok1 + "\t" + tok2, self.alleleSeqs[1:])
        return self.chrom + "\t" + str(self.refPos) + "\t" + self.name + "\t" + self.alleleSeqs[0] + "\t" + altAlleleSeqsStr

