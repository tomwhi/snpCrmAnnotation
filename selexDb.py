#!/usr/bin/env python

import datetime, numpy, pdb, MySQLdb, os, sys, traceback, utility
import illumina, motifModule


# NOTE: I have written "MARKER" in locations where the code must be changed
# if/when the guide file is changed.


def getTF_Name(ensID, db):
    """Returns the HGNC name corresponding to the specified ensemble ID."""
    cmdStr = "SELECT HGNC_Name from TranscriptionFactors where EnsembleID = \"" + ensID + "\";"
    cursor=db.cursor()
    cursor.execute(cmdStr)
    rows = cursor.fetchall()
    if len(rows) > 0:
        return rows[0][0]
    else:
        return None


def getTF_Families_exp(barcode, batch, db):
    """Returns the Vaquerizas family name corresponding to the TF specified
    by experiment."""
    cmdStr = "SELECT vi.TF_Family from Vaquerizas_InterproAnnot vi inner join Vaquerizas_AllPossible_DBD_InterproIDs va on vi.InterproID = va.InterproID inner join Vaquerizas_TF_DBDs vd on va.InterproID = vd.InterproID inner join Vaquerizas_TFs vt on vd.EnsembleID = vt.EnsembleID inner join TranscriptionFactors tf on vt.EnsembleID = tf.EnsembleID inner join CloneDetail cd on tf.EnsembleID = cd.EnsembleID inner join ProteinProduction pp on cd.CloneID = pp.CloneID inner join SELEX_experiments se on pp.ProductionID = se.ProductionID where se.Barcode = \"" + barcode + "\" and se.BarcodeBatch = \"" + batch +"\";"
    cursor=db.cursor()
    cursor.execute(cmdStr)
    rows = cursor.fetchall()
    if (len(rows) > 0):
        return map(lambda row: row[0], rows)
    else:
        return None


def getTF_Name_exp(barcode, batch, db):
    """Returns the HGNC name corresponding to the specified experiment."""
    cmdStr = "SELECT HGNC_Name from TranscriptionFactors tf inner join CloneDetail cd on tf.EnsembleID = cd.EnsembleID inner join ProteinProduction pp on cd.CloneID = pp.CloneID inner join SELEX_experiments se on pp.ProductionID = se.ProductionID where se.Barcode = \"" + barcode + "\" and se.BarcodeBatch = \"" + batch + "\";"
    cursor=db.cursor()
    cursor.execute(cmdStr)
    rows = cursor.fetchall()
    if len(rows) > 0:
        return rows[0][0]
    else:
        return None


def getEnsembleID(tfName, db):
    """Returns the ensemble ID corresponding to the specified HGNC name."""
    cmdStr = "SELECT EnsembleID from TranscriptionFactors where HGNC_Name = \"" + tfName + "\";"
    cursor=db.cursor()
    cursor.execute(cmdStr)
    rows = cursor.fetchall()
    if len(rows) > 0:
        return rows[0][0]
    else:
        return None


def getRelatedTFs(db, ensID, seqSim=0.75):
    """Retrieves the EnsembleIDs of closely-related TFs for a TF specified by
    EnsembleID"""

    # Query the database to get all (ensID2, similarity) tuples for TF1=ensID
    # where the similarity exceeds the specified threshold...
    
    cmdStr = "SELECT TF2, similarity from TF_Similarity where TF1 = \"" + ensID + "\" and similarity > " + str(seqSim) + ";"
    cursor=db.cursor()
    cursor.execute(cmdStr)
    rows = cursor.fetchall()

    if len(rows) == 0:
        return []

    # Sort according to similarity, with the most similar ones at the start:
    resultTups = list(rows)
    resultTups.sort(simCmp)
    assert (resultTups[0][1] == 1)
    return resultTups[1:] # Use 1: to exclude the TF1 itself


def simCmp(tup1, tup2):
    if (tup1[1] > tup2[1]):
        return -1
    elif (tup1[1] == tup2[1]):
        return 0
    else:
        return 1


def getBestMotifRelated(ensembleID, selDb, tomDb, consCountThresh=10000,
                        seqSim=0.75):
    """Returns the best motif for the specified TF. If no motif exists for it
    in either the selex or jaspar databases, then returns the best motif of
    the most closely-related TF that has a motif. If none of the closely
    related TFs have a motif then it returns None."""

    bestMotif = getBestMotif(ensembleID, selDb, tomDb)

    if bestMotif != None:
        return bestMotif
    else:
        # Get the ensembleIDs of all closely-related TFs, ordered according to
        # similarity:
        ensIDs_sims = getRelatedTFs(selDb, ensembleID, seqSim=seqSim)

        # While no motif has been found and there is another closely-related TF
        # not yet considered...
        tfIdx = 0
        while (tfIdx < len(ensIDs_sims)):
            currTup = ensIDs_sims[tfIdx]
            currRelEnsID = currTup[0]
            currRelSim = currTup[1]

            # Retrieve the best motif from the current most closely-related TF:
            currTFbestMotif = getBestMotif(currRelEnsID, selDb, tomDb)
            if (currTFbestMotif != None):
                # The current closely-related TF has a motif -> Return it:
                return currTFbestMotif
            tfIdx = tfIdx + 1
    return None


def getBestMotif(ensembleID, selDb, tomDb, consCountThresh=10000):
    """Returns a single best motif for the TF with the specified ensembleID."""

    (selexMotifs, jasparMotifs) = getMotifs(ensembleID, selDb, tomDb, consCountThresh=consCountThresh)

    if len(selexMotifs) > 0:
        # Consider each selex motif retrieved...
        currMotifIdx = 0
        bestMotif = selexMotifs[currMotifIdx]
        bestMotifConsCount = bestMotif.getConsCount()
        while (currMotifIdx < len(selexMotifs)):
            if (selexMotifs[currMotifIdx].getConsCount() > bestMotifConsCount):
                bestMotif = selexMotifs[currMotifIdx]
                bestMotifConsCount = bestMotif.getConsCount()
            currMotifIdx = currMotifIdx + 1
        return bestMotif
    elif len(jasparMotifs) > 0:
        # No selex motifs were retrieved. Consider jaspar motifs instead...

        # Just arbitrarily return the last jaspar motif found:
        return jasparMotifs[-1]
    else:
        return None


def getAllMotifsRelated(ensID, selDb, tomDb, seqSim=0.75, cycle=3, width=10):
    """Returns a list of pwm objects corresponding to the motifs stored for the
    transcription factor with the specified ensemble ID."""

    # Get ensembleIDs for all TFs to be considered:
    relatedTFs = map(lambda tup: tup[0], getRelatedTFs(selDb, ensID, seqSim=seqSim))
    relatedTFs.append(ensID)

    outputMotifs = []

    for tfID in relatedTFs:
        # Get all the motifs for the current TF:
        currTF_motifs = getMotifs(tfID, selDb, tomDb, consCountThresh=0,
                                  cycle=cycle, width=width)
        currTF_motifsList = currTF_motifs[0]
        currTF_motifsList = currTF_motifsList + currTF_motifs[1]
        outputMotifs = outputMotifs + currTF_motifsList
    return outputMotifs


def getMotifs(ensID, selDb, tomDb, consCountThresh=10000, cycle=3, width=10,
              filteredCons="CCCACCC"):
    """Returns a tuple containing two lists; all SELEX motifs and all JASPAR
    motifs, associated with the specified ensembleID. The selex motifs are
    objects of the type SELEX_pwm. The jaspar motifs are objects of the type
    JASPAR_pwm."""

    # Retrieve all SELEX pwms for this ensembleID...

    # Construct and run a mysql query to retrieve all relevant motifIDs from the
    # SELEX.Motifs table:
    cmdStr = "SELECT distinct m.MotifID from Motifs m inner join SELEX_Analyses sa on m.AnalysisID = sa.AnalysisID inner join Quality q on sa.BarcodeBatchAnalysed = q.Batch and sa.BarcodeAnalysed = q.Barcode where m.EnsembleID = \"" + ensID + "\" and m.ConsensusCount > " + str(consCountThresh) + " and q.Quality = \"ok\" and m.SelexCycle = " + str(cycle) + " and m.width = " + str(width) + ";"
    cursor=selDb.cursor()
    cursor.execute(cmdStr)
    rows = cursor.fetchall()
    motifIDs = map(lambda tup:tup[0], rows)
    
    selexMotifs = []

    # For each motif ID found...
    for motifID in motifIDs:
        # Generate a new motif from this motif's data in the database:
        currMotif = motifModule.SELEX_pwm(selDb, motifID)

        # Filter out those with unwanted consensus:
        if (currMotif.getConsensus() != filteredCons):
            selexMotifs.append(currMotif)

    # Retrieve all JASPAR pwms for this ensembleID...

    # Get all IDs of motifs with this ensembleID:
    cmdStr = "select m.ID from MATRIX m inner join kgAlias2ensembleID kga on m.NAME=kga.kgAlias where kga.ensID = \"" + ensID + "\";"
    cursor=tomDb.cursor()
    cursor.execute(cmdStr)
    rows = cursor.fetchall()
    jasparIDs = map(lambda tup:tup[0], rows)

    jasparMotifs = []

    # For each of those jaspar motifs found...
    for jasparID in jasparIDs:
        # Generate a new motif from this motif's data in the database:
        currMotif = motifModule.JASPAR_pwm(tomDb, jasparID)
        jasparMotifs.append(currMotif)

    return (selexMotifs, jasparMotifs)


def getExpressedTFs(db, tissue, specific=True, minExpn=200, specFoldDiff=10):
    """Returns a list of ensembleIDs for all TFs that are expressed in
    the specified tissue. If "specific", then it filters for only those TFs
    that are *specifically* expressed, using the specified fold difference
    wrt the TF's median expression level across all tissues."""

    if not specific:
        # Run a mysql query to get all TF ensembleIDs where expnLevel > minExp:
        cmdStr = "SELECT distinct tf.EnsembleID from TranscriptionFactors tf inner join geneProbes gp on (tf.EnsembleID = gp.EnsembleID) inner join GNF gnf on gp.probeset = gnf.probeset inner join Vaquerizas_TFs vt on tf.EnsembleID = vt.EnsembleID where gnf.tissue = \"" + tissue + "\" and gnf.expnLevel > " + str(minExpn) + ";"
        cursor=db.cursor()
        print >> sys.stderr, cmdStr
        cursor.execute(cmdStr)
        rows = cursor.fetchall()
        allExpressedTFs = map(lambda tup:tup[0], rows)
        return allExpressedTFs
    else:
        allSpecExpressedTFs = []

        # Run a mysql query to get all TF probeset, expnLevel pairs where
        # expnLevel > minExp:
        cmdStr = "SELECT distinct gp.probeset, gnf.expnLevel from TranscriptionFactors tf inner join geneProbes gp on (tf.EnsembleID = gp.EnsembleID) inner join GNF gnf on gp.probeset = gnf.probeset inner join Vaquerizas_TFs vt on tf.EnsembleID = vt.EnsembleID where gnf.tissue = \"" + tissue + "\" and gnf.expnLevel > " + str(minExpn) + ";"
        cursor=db.cursor()
        cursor.execute(cmdStr)
        allProbesetsExpns = cursor.fetchall()

        # Filter for those TFs whose expression level in this tissue is
        # specFoldDiff times their median expression level over all tissues...
        for (probeset, expn) in allProbesetsExpns:
            cmdStr = "SELECT expnLevel from GNF where probeset = \"" + probeset + "\";"
            cursor=db.cursor()
            cursor.execute(cmdStr)
            rows = cursor.fetchall()
            allExpnLevels = map(lambda tup:tup[0], rows)
            medianExpnLevel = numpy.median(allExpnLevels)
            if (expn > specFoldDiff*medianExpnLevel):
                # This tissue's expression level exceeds the specified fold
                # difference times the median tissue expression level =>
                # this probeset is "specifically" expressed in this tissue.
                # Obtain the probeset's ensemble IDs and add them to the list
                # of specifically expressed TFs:
                cmdStr = "SELECT distinct EnsembleID from geneProbes where probeSet = \"" + probeset + "\";"
                cursor=db.cursor()
                cursor.execute(cmdStr)
                rows = cursor.fetchall()
                currProbesetEnsIDs = map(lambda tup:tup[0], rows)
                allSpecExpressedTFs = allSpecExpressedTFs + currProbesetEnsIDs
            # Uniquify the ensemble IDs of specifically-expressed TFs:
            specExpressedDict = {}
            for tf in allSpecExpressedTFs:
                if (not specExpressedDict.has_key(tf)):
                    specExpressedDict[tf] = 1
        return specExpressedDict.keys()


def loadAnalysisIntoDb(batch, analysisType, command, date, comment, outputDir,
                       barcode, db):
    """Loads specified data into the SELEX_Analyses table."""
    
    # Generate a mysql command for loading the specified parameters into the
    # database:
    cmdStr = "INSERT INTO SELEX_Analyses (BarcodeBatchAnalysed, AnalysisType, CommandUsed, DateRun, Comment, OutputLocation, BarcodeAnalysed) VALUES (\"" + batch + "\",\"" + analysisType + "\",\"" + command + "\",\"" + date + "\",\"" + comment  + "\",\"" + outputDir + "\",\"" + barcode + "\");"
    
    # Run that command:
    cursor=db.cursor()
    cursor.execute(cmdStr)

    # Obtain the new autoincremented analysisID value:
    newAnalysisID = db.insert_id()
    return newAnalysisID


def experimentExists(batch, barcode, db):
    """Returns true if there is such an experiment in the SELEX_experiments
    table, false otherwise."""
    cmdStr = "SELECT * from SELEX_experiments where BarcodeBatch = \"" + \
            batch + "\" and Barcode = \"" + barcode + "\";"
    cursor=db.cursor()
    cursor.execute(cmdStr)
    rows = cursor.fetchall()
    if (len(rows) > 0):
        return True
    else:
        return False


class SelexSample:
    """Encapsulates a single cycle from a single selex experiment. NOTE:
    This class is similar to the older SELEX_sampleInfo class below. However,
    I decided to keep them as separate classes rather than merging into one."""

    def __init__(self, barcode, batch, cycle, passwd="", db="SELEX",
                 user=""):
        self.barcode = barcode
        self.batch = batch
        self.cycle = cycle
        # Database that contains other information about this sample:
        self.db=MySQLdb.connect(passwd=passwd,db=db,user=user)

    def getNpos(self, kLen):
        """Returns the number of positions for sequences of length kLen in
        the sequences for this selex sample."""
        cursor=self.db.cursor()
        cmdStr = "SELECT Total from KmerCountsTotal where Barcode " + \
            " = \"" + self.barcode + "\" and BarcodeBatch = \"" + \
            self.batch + "\" and CycleNo = " + self.cycle + \
            " and K = " + str(kLen) + ";"
        cursor.execute(cmdStr)
        rows = cursor.fetchall()
        return rows[0][0]

    def getHugoName(self):
        """Returns the HUGO gene name for the transcription factor studied by
        this selex experiment."""

        cursor=self.db.cursor()

        # Formulate a mysql query string for retrieving the ensembleID for this
        # selex experiment, using the barcode and batch attributes:
        cmdStr = "select tf.HGNC_Name from TranscriptionFactors tf inner join CloneDetail c on (tf.EnsembleID = c.EnsembleID) inner join ProteinProduction p on (p.CloneID = c.CloneID) inner join SELEX_experiments se on (se.ProductionID = p.ProductionID) where (se.Barcode = \"" + self.barcode + "\" and se.BarcodeBatch = \"" + self.batch + "\");"
        cursor.execute(cmdStr)
        rows = cursor.fetchall()

        # Return the first element of the first result from the query run:
        if (len(rows) > 0):
            return rows[0][0]
        else:
            # No hugo gene name was available:
            return None

    def getEnsembleID(self):
        """Returns the ensembleID for the transcription factor studied by
        this selex experiment."""

        cursor=self.db.cursor()

        # Formulate a mysql query string for retrieving the ensembleID for this
        # selex experiment, using the barcode and batch attributes:
        cmdStr = "select c.EnsembleID from CloneDetail c inner join ProteinProduction p on (p.CloneID = c.CloneID) inner join SELEX_experiments se on (se.ProductionID = p.ProductionID) where (se.Barcode = \"" + self.barcode + "\" and se.BarcodeBatch = \"" + self.batch + "\");"
        cursor.execute(cmdStr)
        rows = cursor.fetchall()

        # Return the first element of the first result from the query run:
        if (len(rows) > 0):
            return rows[0][0]
        else:
            # No ensembleID was available:
            return None

    def getName(self):
        """Returns a string representation of this selex sample."""

        return self.barcode + "_" + self.batch + "_" + \
            str(self.cycle) + "_" + str(self.getHugoName())

    def getSampleID(self):
        """Returns a string representation of this selex sample's ID."""

        return self.barcode + "_" + self.batch + "_" + \
            str(self.cycle)


class SELEX_sampleInfo:
    """This class contains metadata for a single SELEX sample; i.e.
    a single well in a 96-well plate in which a SELEX experiment was
    conducted."""

    def __init__(self, tokens, parentSeqRunInfo):
        # tokens is a list of string tokens from which this sampleInfo
        # will be generated.
        # parentSeqRunInfo is an object representing the sequencing run
        # that this sample was associated with.

        self.parentSeqRunInfo = parentSeqRunInfo

        # Parse the tokens of the guide file into the information that
        # will ultimately be input into the MySQL database...
        # MARKER: This section must change if the ordering/contents of elements
        # in the guide file is modified...
        plateAndWell = tokens[1].split("_")
        self.plate = int(plateAndWell[0])
        self.well = plateAndWell[1]
        self.barcodeName = tokens[2]
        self.regExp = tokens[3]
        self.substitute1 = tokens[4]
        self.substitute2 = tokens[5]
        self.batch = tokens[6]
        self.arttusGeneName = tokens[7]
        self.cycle = int(tokens[8])

        # NOTE: The ensemble ID must match the pattern "ENS.*", and
        # must be for a factor already in the transcription factor
        # annotation database, otherwise it will not link to the
        # transcription factor annotation tables.
#            self.ensID = None
#            if (tokens[9][:3] == "ENS"):
        self.ensID = tokens[9]
        self.domain = tokens[10]
        self.type = tokens[11]
        self.namedOrganism = tokens[12]
        self.production = tokens[13]
        self.expression = tokens[14]
        self.note = tokens[15]

        # tokens at indeces 16 and 17 are for flowcell and lane,
        # respectively.
        self.aaSeq = tokens[18]
        self.sequenced = tokens[19]

        # Check that global tokens (that must be the same for all samples for
        # this sequencing run - e.g. FlowcellID, illuminaFilename etc.) are
        # identical to the corresponding fields in parentSeqRunInfo:
        # MARKER: This bit potentially needs to be updated too:
        if (tokens[0] != self.parentSeqRunInfo.illuminaFilename):
            raise ValueError("Illumina file name is internally inconsistent.")
        if (tokens[-4] != self.parentSeqRunInfo.flowcell):
            raise ValueError("Specified flow cell is internally inconsistent.")
        if (tokens[-3] != self.parentSeqRunInfo.lane):
            raise ValueError("Specified lane is internally inconsistent.")


class SELEX_seqRunInfo:
    """This class contains metadata for a single sequencing run for a set of
    SELEX analyses."""

    def __init__(self, guideFileName, qualThresh=20):
        # The name of the tab-delimited file containing the metadata:
        self.guideFileName = guideFileName

        # Retrieve the "global" information (pertaining to all samples
        # associated with this sequencing run) from the file...
        guideFile = open(self.guideFileName)

        # The global metadata is the illumina filename, the flowcellID, and
        # the lane. MARKER: If the guide file structure changes, then this
        # section may have to change too:
        firstDataLine = guideFile.readlines()[1]
        toks = firstDataLine.strip().split("\t")

        self.illuminaFilename = toks[0]
        self.flowcell = toks[-4]
        self.lane = toks[-3]

        # NOTE: Will record the date when the info is actually loaded up, rather
        # than recording it now.

        # Parse the metadata for the individual samples:
        self.sampleInfos = self.parseSampleInfos()

        # MOVED THIS FROM checkSeqData() method to here, on April 6th 2011:
        # Make a new oligoDesignFrequency object, from the list of oligoDesign
        # regular expressions in self.sampleInfos...
        # NOTE: Modified on May 27th 2011: Storing barcode name information
        # too now...
        designObjs = []
        for sampleInfo in self.sampleInfos:
            designStr = sampleInfo.regExp
            barcodeName = sampleInfo.barcodeName

            # Extract the barcode, nRandom, and constEndStr from the current
            # oligo design regExp, and make a new oligoDesign object
            # encapsulating this information:
            barcode = designStr[1:designStr.find(".")]
            endSection = designStr[designStr.rfind(".")+1:]
            numberNs = \
                len(designStr[designStr.find("."):designStr.rfind(".")+1])
            currDesign = oligoDesign(None, barcodeName, barcode, numberNs,
                                     endSection)
            designObjs.append(currDesign)
        self.designSet = oligoDesignSet(designObjs, qualThresh=qualThresh)

    def getDesignSet(self):
        return self.designSet

    def parseSampleInfos(self):
        sampleInfos = []
        lineNum = 2
        for line in open(self.guideFileName).readlines()[1:]:
            elems = line.strip("\n").split("\t")

            # Ignore empty lines:
            if (len(elems) > 1):
                try:
                    sampleInfos.append(SELEX_sampleInfo(elems, self))
                except ValueError, e:
                    print >> sys.stderr, \
                        "Invalid SELEX sample information; ValueError " + \
                        " occurred, line", lineNum, "of input file:"
                    print >> sys.stderr, line
                    print >> sys.stderr, "ValueError text:"
                    print >> sys.stderr, e
                    sys.exit(1)
                except IndexError, e:
                    print >> sys.stderr, \
                        "Invalid SELEX sample information; IndexError " + \
                        " occurred, line", lineNum, "of input file:"
                    print >> sys.stderr, line
                    print >> sys.stderr, "IndexError text:"
                    print >> sys.stderr, e
                    sys.exit(1)
            lineNum = lineNum + 1

        return sampleInfos

    def getSampleInfos(self):
        return self.sampleInfos

    def checkSeqData(self, outFile, sanger=False):
        """This method does barcode counting and sequence quality counting,
        and reports the outcomes to the specified output file."""

        # Check that the illumina file exists:
        if (not os.path.exists(self.illuminaFilename)):
            print >> sys.stderr, "Specified illumina file \"" + \
                self.illuminaFilename + "\" does not exist."
            sys.exit(1)

        # Count the instances of the designs, and also sequence quality, in
        # the illumina input file:
        print >> sys.stderr, "Counting oligo design frequencies in input" + \
            " illumina dataset " + self.illuminaFilename, "..."
        self.designSet.countOligoDesigns(open(self.illuminaFilename),
                                         sanger=sanger)
        print >> sys.stderr, "Finished counting oligo design frequencies."

        # Print the relevant information in designSet out to the specified
        # output file (i.e. print each design along with it's count, and also
        # print the total number of zero-mapping and multiple-mapping sequence
        # reads)...

        # I will print the information in the following format...

        # First, the total number of sequences:
        print >> outFile, "numSeqs = ", self.designSet.nSeqsTotal

        # Then, the number of sequences of poor quality (along with the specified quality threshold used for determining this):
        print >> outFile, "numPoorQuality = ", self.designSet.nPoorQualitySeqs

        # Then, the fraction (calculated from the above two values):
        print >> outFile, "fracPoorQuality = ", self.designSet.nPoorQualitySeqs/float(self.designSet.nSeqsTotal)

        print >> outFile, "numGoodQuality = ", \
            self.designSet.nSeqsTotal - self.designSet.nPoorQualitySeqs

        # Then, the number of zero-"mapping" sequences:
        print >> outFile, "numZeroMappers = ", self.designSet.numNonMatching

        # Then, the number of multi-"mapping" sequences:
        print >> outFile, "numMultiMappers = ", self.designSet.numMultipleMatching

        # Then, for each oligo design, print the design (barcode, gap, end) and
        # number of sequences matching that design...
        for design in self.designSet.getDesigns():
            print >> outFile, "Count of design", design.toString(), \
                "=", design.getCount()


class oligoDesign:
    """This class encapsulates the concept of an oligo design. It is used
    to facilitate counting of the number of occurrences of each design in
    a sequence file (at least that is what I intend it for)."""

    def __init__(self, designID, barcodeName, barcode, nRandom, constEndStr):
        self.designID = designID
        # NOTE: Added this attribute on May 27th 2011. This represents the
        # actual barcode string from the guidefile, as opposed to the barcode
        # value extracted from the regular expression field in the same file:
        self.barcodeName = barcodeName
        self.barcode = barcode
        self.nRandom = nRandom
        self.constEndStr = constEndStr

        # Optional file that will contain sequences matching this design:
        self.seqsFile = None

        # Optional number of occurrences of this design in some given set of
        # sequences:
        self.count = 0

        # Optional string describing a factor that this oligo design was
        # used to study in a particular sequencing run:
        self.factorString = None

    def getBarcodeName(self):
        return self.barcodeName

    def getID(self):
        return self.designID

    def incrementCount(self):
        self.count = self.count + 1

    def setFactorString(self, factorString):
        self.factorString = factorString

    def getSeqsFile(self):
        return self.seqsFile

    def getRanLength(self):
        return self.nRandom

    def getRandomSection(self, seq):
        """Extracts the random section from the specified sequence, assuming
        it exhibits this oligo design. Returns that string."""
        return seq[:(len(self.getBarcode()) + self.getRanLength())]\
            [-self.getRanLength():]

    def getCount(self):
        return self.count

    def setSeqsFile(self, seqsFile):
        self.seqsFile = seqsFile

    def getBarcode(self):
        return self.barcode

    def matches(self, sequence):
        """Reports whether or not the specified sequence matches this design."""

        # Check barcode first:
        if (sequence[:len(self.barcode)] != self.barcode):
            return False

        # Check constant end string:
        constStrStart = len(self.barcode) + self.nRandom
        constStrEnd = constStrStart + len(self.constEndStr)
        sequenceFlankSection = sequence[constStrStart:constStrEnd]
        if (sequenceFlankSection != self.constEndStr):
            return False

        # The sequences matches this oligo design:
        return True

    def toString(self):
        """Returns a string representation of this object."""
        return self.barcode + "(" + str(self.nRandom) + "N)" + self.constEndStr


class oligoDesignSet:
    """This class represents a set of oligo designs. This facilitates such
    things as recording the counts of oligo designs expected to be in a
    given illumina file, and also extracting sequence data for a set of
    oligo designs that were analysed in a sequencing run."""

    def __init__(self, oligoDesigns, qualThresh=20):
        # oligoDesigns is a list of oligoDesign objects

        # Determine all unique length values in the specified designs:
        lengths = []
        for design in oligoDesigns:
            if design.getRanLength() not in lengths:
                lengths.append(design.getRanLength())

        # List of unique random section lengths in the set of designs:
        self.lengths = lengths

        # Determine the minimum barcode length amongst all oligo designs
        # specified:
        self.minBarcodeLength = min(map(lambda design: len(design.barcode),
                                        oligoDesigns))

        # Set up a hashtable to facilitate quick mapping of an individual
        # sequence to a single oligo design or set of oligo designs.
        # Each key is a minimum length barcode string (minimum length because
        # different barcodes can be different lengths), and each value is a
        # list of all designs whose barcode start sequence is the same as that
        # minimum length barcode string...
        self.designs = {}
        for currOligoDesign in oligoDesigns:
            barcodeKey = currOligoDesign.barcode[:self.minBarcodeLength]
            if (self.designs.has_key(barcodeKey)):
                self.designs[barcodeKey].append(currOligoDesign)
            else:
                self.designs[barcodeKey] = [currOligoDesign]

        # The following attributes record optional information about
        # the occurrence of some set of sequences wrt this set of
        # oligo designs...

        # Phred score cutoff for single worst-quality base from a sequence:
        self.seqQualThresh = qualThresh

        # Number of sequences matching no barcode:
        self.numNonMatching = 0
        # Number of sequences matching multiple barcodes:
        self.numMultipleMatching = 0

        # Will store info about sequence quality count, even though it
        # doesn't fit in very well conceptually... If I really wanted
        # to make this conceptually clean, perhaps I would use a subclass,
        # but presently that seems like overkill.
        self.nPoorQualitySeqs = 0
        # Will also count the total number of sequences in the file:
        self.nSeqsTotal = 0

    def getLengths(self):
        """Returns an array containing all the unique random-section lengths
        that are observed in all of the oligo designs in this set."""
        return self.lengths

    def db2SeqFiles(self, db, seqRunID, randomOnly=True):
        """Writes all the sequences corresponding to the specified sequencing
        run ID (in the specified selex database) to individual output files for
        each oligoDesign in this set. Reports (to stderr) number of zero-mapping
        and multi-mapping sequences observed."""

        # Retrieve all the specified sequences and write them to a temporary
        # sequence file...
        
        # Place the temporary file in "/tmp":
        # Modified on June 8th 2011: Writing to /home/tom/tom/mysqlTmpFiles
        # instead of /tmp. This allows me to delete the temporary file after
        # I'm finished with it:
        tmpSeqsFilename = \
            utility.makeTempFilename("/home/tom/tom/mysqlTmpFiles/db2SeqFiles_" + str(seqRunID) + "_")
        writeSeqsToFile(db, seqRunID, tmpSeqsFilename)

        # Output the sequences in that temorary file to the output files
        # referenced by each of the oligo designs in this set:
        self.seqFile2SeqFiles(tmpSeqsFilename, randomOnly=randomOnly)
        utility.deleteFiles([tmpSeqsFilename])

    def illumFile2dbFile(self, illumFilename, tableFile, seqRunID):
        """This method takes as input the name of a file containing illumina
        sequences, and an open file handle to which the mysql-formatted table
        data should be written. It writes out all of the sequences to that table
        file, inferring a barcode string for each sequence and writing that
        information out as well as the other info such as the sequence itself
        and it's quality score."""

        # Create an illumSeqParser object using the illumina file as input:
        seqParser = illumina.illumSeqParser(open(illumFilename))

        # For each sequence in the specified text file of DNA sequences...
        seqCount = 0
        while seqParser.hasSeq():
            # Update progress message:
            if ((seqCount % 1000000) == 0):
                print >> sys.stderr, "Processed", seqCount, "sequences."
            seqCount = seqCount + 1
            currSeq = seqParser.nextSeq()

            # Try to retrieve the oligoDesign that the sequence exhibits:
            matchingDesigns = self.getMatchingDesigns(currSeq.getSeq())

            barcodeStr = None
            if (len(matchingDesigns) > 1):
                # Multiple designs match. Generate a barcode string that
                # reflects this fact, for this sequence:
                barcodeStr = "MultiMatch"

                # April 14th 2011: Won't record which designs this sequence
                # matched to anymore:
                #for matchingDesign in matchingDesigns:
                #    barcodeStr = barcodeStr + "_" + matchingDesign.toString().replace("(","").replace(")","")
            elif (len(matchingDesigns) == 0):
                # This sequence does not match any of the designs =>
                # Generate a barcode string that reflects this:
                barcodeStr = "NoMatch"
            else:
                # There was exactly one match => Obtain the barcode string
                # for that match:
                barcodeStr = matchingDesigns[0].getBarcodeName()

            # Output this sequence information to the specified output
            # database table file:
            outStr = '"' + barcodeStr + '","' + currSeq.sequence + '","' + \
                currSeq.quality + '",' + str(seqRunID)
            print >> tableFile, outStr

    def seqFile2SeqFiles(self, seqFilename, randomOnly=True):
        """Redirects all sequences in the specified text file to the text files
        of matching oligo designs in this set of oligo designs."""

        # For each sequence in the specified text file of DNA sequences...
        seqCount = 0
        for line in open(seqFilename).xreadlines():
            # Update progress message:
            if ((seqCount % 1000000) == 0):
                print >> sys.stderr, "Processed", seqCount, "sequences."
            seqCount = seqCount + 1

            currSeq = line.strip()

            # Try to retrieve the oligoDesign that the sequence exhibits:
            matchingDesigns = self.getMatchingDesigns(currSeq)

            if (len(matchingDesigns) > 1):
                self.numMultipleMatching = self.numMultipleMatching + 1
            elif (len(matchingDesigns) == 0):
                self.numNonMatching = self.numNonMatching + 1
            else:
                # There was exactly one match => Write the sequence to the
                # single mapped oligo design's file...

                # If only retrieving the random section (i.e. randomOnly==True),
                # then get the random section from the matching design:
                if (randomOnly):
                    outSeq = matchingDesigns[0].getRandomSection(currSeq)
                else:
                    # Otherwise will output the entire sequence:
                    outSeq = currSeq

                # Write the sequence out:
                matchingDesigns[0].getSeqsFile().write(outSeq + "\n")
                matchingDesigns[0].getSeqsFile().flush()

        # Make sure all of the output above is written out without delay...
        for currDesign in self.getDesigns():
            currDesign.getSeqsFile().flush()

    def getDesigns(self):
        """Returns a list of all oligoDesign objects in this set."""
        allDesigns = []
        for designArr in self.designs.values():
            allDesigns = allDesigns + designArr
        return allDesigns

    def getMatchingDesigns(self, sequence):
        """Determines which oligo designs in this set match the specified DNA
        sequence. Returns an array of those designs."""

        # Retrieve the list of all oligo designs that might match the sequence,
        # according to the starting section of the sequence (minimum barcode
        # length in size):

        seqMinBarcodeSection = sequence[:self.minBarcodeLength]
        if (not self.designs.has_key(seqMinBarcodeSection)):
            return []

        possibleMatches = self.designs[seqMinBarcodeSection]
        designMatches = []
        for design in possibleMatches:
            if design.matches(sequence):
                designMatches.append(design)
        return designMatches

    def countOligoDesigns(self, illuminaSeqFile, sanger=False):
        """Counts the frequency of each of the oligo designs, plus the other
        misc sequences and sequences that map to more than one design, in
        sequences from the specified input illumina sequence file.
        Also counts the total number of sequences and the number with
        a quality score below the specified threshold."""

        # NOTE: The sequence quality counting doesn't really fit in here
        # conceptually, but it is done here anyway to avoid having to
        # read and parse the sequence data twice (which is relatively
        # time consuming).

        # Create an illumina file parser to process the input file:
        seqParser = illumina.illumSeqParser(illuminaSeqFile, sanger=sanger)

        seqCount = 1
        while seqParser.hasSeq():
            currSeq = seqParser.nextSeq()
            if ((seqCount % 100000) == 0):
                print >> sys.stderr, "Counted", seqCount, "sequences..."
            seqCount = seqCount + 1

            # Only count barcode matches for those sequences that pass
            # the quality cutoff...
            if (currSeq.getMinQual() >= self.seqQualThresh):
                # Get all matching designs for the current sequence:
                matchingDesigns = self.getMatchingDesigns(currSeq.getSeq())

                # Record the occurrence of zero-mappers and multi-mappers:
                if (len(matchingDesigns) == 0):
                    # No matches were found. Record this fact:
                    self.numNonMatching = self.numNonMatching + 1
                elif (len(matchingDesigns) > 1):
                    # Multiple oligoDesigns were found that match this
                    # sequence. Record this fact:
                    self.numMultipleMatching = self.numMultipleMatching + 1

                # Record the fact that each individual design was observed:
                for design in matchingDesigns:
                    design.incrementCount()

            # The sequence was poor quality. Record this:
            else:
                self.nPoorQualitySeqs = self.nPoorQualitySeqs + 1
            self.nSeqsTotal = self.nSeqsTotal + 1

                # # Get the first minBarcodeSize letters from the sequence, to
                # # faciliate retrieval of oligoDesign objects potentially
                # # matching this sequence:
                # currSeqBarcodeKey = currSeq.getSeq()[:self.minBarcodeLength]

                # # Retrieve the array of oligo designs (and their counts) whose
                # # barcode starts with that sequence. Do this by simply using the
                # # hashtable of this object:
                # if (self.designs.has_key(currSeqBarcodeKey)):
                #     possibleDesignMatches = \
                #         self.designs[currSeqBarcodeKey]
                # else:
                #     possibleDesignMatches = []
                
                # # Consider each possible matching design...
                # nDesignMatches = 0
                # for designCount in possibleDesignMatches:
                #     if designCount[0].matches(currSeq.getSeq()):
                #         # This design matched the sequence; record that fact
                #         designCount[1] = designCount[1] + 1
                #         # Record that an additional matching design was
                #         # found:
                #         nDesignMatches = nDesignMatches + 1
                # if (nDesignMatches == 0):
                #     # No matches were found. Record this fact:
                #     self.numNonMatching = self.numNonMatching + 1
                # elif (nDesignMatches > 1):
                #     # Multiple oligoDesigns were found that match this
                #     # sequence. Record this fact:
                #     self.numMultipleMatching = self.numMultipleMatching + 1


class dataLoader:
    """This class is used to map and load data from a SELEX_seqRunInfo object
    into the local MySQL database.

    NOTE: This doesn't really encapsulate a clear concept and so is not great
    OO code. At some point I may need to re-engineer this module wrt this class.
"""

    def __init__(self, seqRunInfo):
        self.seqRunInfo = seqRunInfo

        # Check that the referenced illumina file and the original guide file
        # still exist...
        if (not os.path.exists(seqRunInfo.guideFileName)):
            print >> sys.stderr, "SELEX metadata guide file \"" + \
                seqRunInfo.guideFileName + "\" no longer exists. Aborting."
            sys.exit(1)
        if (not os.path.exists(seqRunInfo.illuminaFilename)):
            print >> sys.stderr, "SELEX illumina sequencing file \"" + \
                seqRunInfo.illuminaFilename + "\" no longer exists. Aborting."
            sys.exit(1)

        # Open a connection to dataLoader's specified server and database:
        self.db=MySQLdb.connect(passwd="",db="SELEX",user="")
        self.cursor=self.db.cursor()

        # Check that no sequenceing run has been loaded where the guide file
        # or illumina file is the same as that of this sequence run info
        # object...

        cmdStr = "SELECT * from SequencingRuns where GuideFileLocation = \"" + \
            self.seqRunInfo.guideFileName + "\";"
        self.cursor.execute(cmdStr)
        rows = self.cursor.fetchall()
        if (len(rows) > 0):
            print >> sys.stderr, "ERROR: Database already contains row", \
                "for Guide File", self.seqRunInfo.guideFileName
            sys.exit(1)

        cmdStr = "SELECT * from SequencingRuns where IlluminaFileLocation" + \
            " = \"" + self.seqRunInfo.illuminaFilename + "\";"
        self.cursor.execute(cmdStr)
        rows = self.cursor.fetchall()
        if (len(rows) > 0):
            print >> sys.stderr, "ERROR: Database already contains row", \
                "for Illumina File", self.seqRunInfo.illuminaFilename
            sys.exit(1)

    def mapAndLoad(self):
        """This method does three things:
- Firstly, maps the metadata from the SELEX_seqRunInfo object to fields in the database schema, and also fills in other fields of modified tables (e.g. foreign key linking IDs).
- Secondly, creates a temporary file containing the sequence data such that it can be loaded to the MySQL database using an efficient "load data infile" command.
- Thirdly, loads the metadata into the database, after the mapping has been done. NOTE: This step and the first step actually have to occur simultaneously, due to the need for correctly assigning foreign key values after updates occur.
- NOTE: As it is *critical* that the order of the loading is done correctly, I implemet these three steps as one monolithic block of code in this method, to avoid the potential for accidentally incrementing a foreign key ID field *after* linking.
"""
        # MARKER: This method must change if the ordering/contents of elements
        # in the guide file *or* the MySQL schema is modified...

        # Turn off MySQL's autocommit:
        # Start a new MySQL transaction:
        # NOTE: I don't think I need to do these things here, as they are the
        # default action when using MySQLdb. => Doing nothing here.

        # Attempt to map and load all the information in the sequencing run
        # info object into the MySQL database...
        try:
            # Map and load information to seqRuns table...
            
            # Calculate today's date, into string variable todaysDate:
            year = datetime.date.today().year
            month = datetime.date.today().month
            day = datetime.date.today().day
            todaysDate = str(year) + "-" + str(month) + "-" + str(day)

            # NOTE: In this code, these commands implicitly state the mapping
            # from the metadata to the MySQL schema...

            # Generate the command for loading the data into the database:
            cmdStr = "INSERT INTO SequencingRuns (GuideFileLocation, FlowCell, Lane, IlluminaFileLocation, Date) VALUES (\"" + self.seqRunInfo.guideFileName + "\",\"" + str(self.seqRunInfo.flowcell) + "\"," + self.seqRunInfo.lane + ",\"" + self.seqRunInfo.illuminaFilename + "\",\"" + todaysDate + "\");"

            # Run the command, to load the data into the mysql database:
            self.cursor.execute(cmdStr)

            # Obtain the new sequenceRunID value resulting from the new,
            # most-recent insertion into the SequencingRuns table:
            newSeqrunID = self.db.insert_id()
            print >> sys.stderr, "New SequencingRunID =", newSeqrunID

            # Map and load the individual "sample" metadata into the
            # correct tables...
            for sampleInfo in self.seqRunInfo.getSampleInfos():
                # Map and load to cloneDetail table...

                # Generate the command for loading the data into the database:
                # FIXME: Map the amino acid sequence here.
                # NOTE: If the ensID isn't already recorded in
                # the TranscriptionFactors table, then don't specify it:
                cmdStr = "SELECT * from TranscriptionFactors where " + \
                    "EnsembleID = \"" + sampleInfo.ensID + "\";"
                self.cursor.execute(cmdStr)
                tf_matches = self.cursor.fetchall()
                cmdStr = None
                # May 7th 2012: Some sampleInfo.domain values have quotes in
                # them, stuffing up the mysql query. Removing those here:
                arttusNameValue = sampleInfo.arttusGeneName.replace("\"", "")
                if (len(tf_matches) < 1):
                    print >> sys.stderr, "WARNING: No TF in database for ", \
                        "ensemble ID", sampleInfo.ensID, ". Marking ID as null."
                    cmdStr = "INSERT INTO CloneDetail (ArttusGeneName, DomainStudied, Length, ArttusNamedOrganism, AminoAcidSequence) VALUES (\"" + arttusNameValue + "\",\"" + sampleInfo.domain + "\",\"" + sampleInfo.type + "\",\"" + sampleInfo.namedOrganism + "\",\"" + sampleInfo.aaSeq +"\");"
                else:
                    cmdStr = "INSERT INTO CloneDetail (EnsembleID, ArttusGeneName, DomainStudied, Length, ArttusNamedOrganism, AminoAcidSequence) VALUES (\"" + sampleInfo.ensID + "\",\"" + arttusNameValue + "\",\"" + sampleInfo.domain + "\",\"" + sampleInfo.type + "\",\"" + sampleInfo.namedOrganism + "\",\"" + sampleInfo.aaSeq + "\");"
                # Run the command, to load the data into the mysql database:
                print >> sys.stderr, cmdStr
                self.cursor.execute(cmdStr)

                # Obtain the resulting cloneID value:
                newCloneID = self.db.insert_id()
                print >> sys.stderr, "New CloneID =", newCloneID
                
                # Map and load to proteinProduction table...
                cmdStr = "INSERT INTO ProteinProduction (ProductionMethod, CloneID, ExpressionLevel) VALUES (\"" + sampleInfo.production + "\"," + str(newCloneID) + ",\"" + sampleInfo.expression + "\");"

                # Run the command, to load the data into the mysql database:
                self.cursor.execute(cmdStr)

                # Obtain the resulting proteinProductionID value:
                newProteinProductionID = self.db.insert_id()
                print >> sys.stderr, "New ProteinProductionID =", \
                    newProteinProductionID

                # Map and load to oligoDesigns table...
                # NOTE: This is possibly where I would decompose the reg-exp
                # strings into more meaningful components, once I understand
                # the full significance of these better.
                cmdStr = "INSERT INTO OligoDesigns (BarcodeName, RegularExpression, Substitute1, Substitute2) VALUES (\"" + sampleInfo.barcodeName + "\",\"" + sampleInfo.regExp + "\",\"" + sampleInfo.substitute1 + "\",\"" + sampleInfo.substitute2 + "\");"

                # Run the command, to load the data into the mysql database:
                self.cursor.execute(cmdStr)
                
                # Obtain the resulting oligoDesignsID value:
                newOligoDesignsID = self.db.insert_id()
                print >> sys.stderr, "New OligoDesignsID =", \
                    newOligoDesignsID

                # Optionally map and load to experiments table, depending on
                # whether an experiment has already been created for this
                # barcode-batch combination...
                # Something else to consider here: Not all (Barcode,Batch)
                # pairs are part of an experiment. In particular, when a sample
                # (with a given (Barcode, Batch) pair) has been marked as
                # "contamination", then it means that sample was not selected
                # for sequencing. Also, if a sample is for cycle zero, then
                # it will not belong to a *single* experiment but will rather
                # be shared by *all* experiments.

                # FIXME: I need to specify the barcode for this table, but i
                # haven't yet inferred the actual barcode itself! At the moment,
                # I will just use sample.barcodeName.

                # Find out whether an experiment with this barcode/batch
                # combination has already been run:
                cmdStr = "SELECT * from SELEX_experiments where " + \
                    "BarcodeBatch = \"" + sampleInfo.batch + "\"" + \
                    "and Barcode = \"" + sampleInfo.barcodeName + "\";"
                self.cursor.execute(cmdStr)
                experiment_matches = self.cursor.fetchall()

                if (len(experiment_matches) < 1):
                    cmdStr = "INSERT INTO SELEX_experiments (BarcodeBatch, Barcode, ProductionID, OligoDesignID) VALUES (\"" + sampleInfo.batch + "\",\"" + sampleInfo.barcodeName + "\"," + str(newProteinProductionID) + "," + str(newOligoDesignsID) + ");"

                    # Run the command, to load the data into the mysql database:
                    self.cursor.execute(cmdStr)

                # Changed on Feb 18th 2011: Here, checking to see whether the
                # sample has already been loaded, and not loading if this is
                # the case. This check becomes necessary once multiple
                # different sequencing lanes are loaded into the database.

                # Changed on Feb 21st 2011: Modified to just load the
                # barcode/batch/cycle info to the "samples" table, and to
                # add the extra info (Note, Plate, Well, SelectedForSequencing)
                # to the SampleSequencings table.

                cmdStr = "SELECT * from SELEX_samples where " + \
                    "BarcodeBatch = \"" + sampleInfo.batch + "\"" + \
                    "and Barcode = \"" + sampleInfo.barcodeName + "\"" + \
                    "and CycleNo = " + str(sampleInfo.cycle) + ";"
                self.cursor.execute(cmdStr)
                sample_matches = self.cursor.fetchall()

                if (len(sample_matches) < 1):
                    # There is no pre-existing sample of this identifier =>
                    # map and load to SELEX_samples table...
                    cmdStr = "INSERT INTO SELEX_samples (BarcodeBatch, Barcode, CycleNo) VALUES (\"" + sampleInfo.batch + "\",\"" + sampleInfo.barcodeName + "\"," + str(sampleInfo.cycle) + ");"

                    # Run the command, to load the data into the mysql database:
                    self.cursor.execute(cmdStr)

                # There is no pre-existing sample of this identifier =>
                # map and load to SELEX_samples table...
                # NOTE: The "note" field sometimes contains quotation
                # characters. This renders the automatically-generated mysql
                # command invalid. Therefore, I replace quotation marks
                # in that string. I should really do this more comprehensively
                # for the other fields, but this would be tricky, given
                # the way I have implemented this method.
                noteStr = sampleInfo.note.replace("\"","\\\"")

                cmdStr = "INSERT INTO SampleSequencings (BarcodeBatch, Barcode, CycleNo, Note, Plate, Well, SequencingRunID, SelectedForSequencing) VALUES (\"" + sampleInfo.batch + "\",\"" + sampleInfo.barcodeName + "\"," + str(sampleInfo.cycle) + ",\"" + noteStr + "\"," + str(sampleInfo.plate) + ",\"" + sampleInfo.well + "\",\"" + str(newSeqrunID) + "\"," + sampleInfo.sequenced + ");"
                
                # Run the command, to load the data into the mysql database:
                self.cursor.execute(cmdStr)

            # Finally, load the actual sequence data into the "Sequences"
            # table...

            # Generate the sequence file in "/tmp", by parsing the sequences
            # file and writing out the resulting sequence data plus the
            # sequenceRunID value obtained from the database...

            # First, check that a file of that name does not already exist:
            guideFilePathEnd = \
                self.seqRunInfo.guideFileName.split("/")[-1].split(".")[0]
            tmpSeqsFilename = "/tmp/" + guideFilePathEnd + "_tmpSeqsTable.txt"
            if (os.path.exists(tmpSeqsFilename)):
                print >> sys.stderr, "Temporary illumina sequences table \"" + \
                    tmpSeqsFilename + "\" already exists. Aborting."
                sys.exit(1)
            
            # INTRODUCED ON APRIL 6TH 2011, TO REPLACE DIRECT USE OF
            # ILLUMPARSER. THIS IS BECAUSE I NOW NEED TO WRITE THE BARCODE
            # STRING TO THE DATABASE AS WELL:
            # Write the sequences out. Use the sequencing run's oligoDesignSet
            # object to achieve this:
            seqsOutfile = open(tmpSeqsFilename, 'w')
            self.seqRunInfo.getDesignSet().illumFile2dbFile(\
                self.seqRunInfo.illuminaFilename, seqsOutfile, newSeqrunID)
            seqsOutfile.flush()
            seqsOutfile.close()

            # Run the "load data infile" mysql command to load the resulting
            # file of data into the sequences table:
            cmdStr = "SET sql_mode='NO_BACKSLASH_ESCAPES'";
            self.cursor.execute(cmdStr)

            # XXX FIXME: Do the load without "local" keyword, since it is
            # much faster, apparently.
            cmdStr = "LOAD DATA LOCAL INFILE '" + tmpSeqsFilename + \
                "' INTO TABLE Sequences FIELDS TERMINATED BY ',' " + \
                "ENCLOSED BY '\"' " + \
                "(Barcode, SeqValue, QualityScore, SequencingRunID);"

            self.cursor.execute(cmdStr)

        except Exception, e:
            # Except all exceptions. If anything went wrong when attempting
            # to load data to the database, then we must "rollback" all changes,
            # report the failure and exit:
            self.db.rollback()
            print >> sys.stderr, "ERROR: Mapping and loading data failed " + \
                "for some reason. Rolling back all changes and quitting."
            print >> sys.stderr, "Printing stack traceback of the exception..."
            print >> sys.stderr, traceback.print_exc()
            sys.exit(1)

        # Commit all the changes that have been made, ending the current
        # transaction...
        self.db.commit()


def writeSeqsToFile(conn, seqRunID, outFileName):
    """Writes all sequences for the specified sequencing run ID to the
    output file with the specified name.

    "conn" is a MysqlDB database connection object for the selex database of
    interest."""

    # Write all the sequences for that sequencing run ID to the specified output
    # file, using a mysql command...
    
    cmdStr = "SELECT SeqValue FROM Sequences WHERE SequencingRunID = " + \
                    str(seqRunID) + " INTO OUTFILE \"" + outFileName \
                    + "\";"
    conn.cursor().execute(cmdStr)

