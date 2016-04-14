#!/usr/bin/env python

import copy, math, numpy, os, pdb, sys
import chipseq_analysis, cmdlineProgs, random, selexDb, seq_utility, seqUtils, utility
import numpy

# SUMMARY: I have created this class to encapsulate motifs, for use in my
# chipseq qc work. Should possibly have done this earlier, as there are
# various motif-related functions etc. scattered throughout my python libraries.


# Introduced on June 18th 2013, for use in snv analysis. Should be useful
# elsewhere too:
def motifPrefixes2meme(outputFile, motifPrefixes, baseDir, suffix):
    for motifPrefix in motifPrefixes:
        currMotif = pwm(open(baseDir + motifPrefix + suffix))
        currMotif.setName(motifPrefix)
        currMotif.writeToMEME(outputFile)
        print >> outputFile, ""


# Introduced on June 18th 2013, for use in snv analysis. Should be useful
# elsewhere too:
def getTFsForMotifPrefix(motifPrefix, motifsCursor):
    """Retrieves all TFs corresponding to the specified motif prefix,
    by identity and sequence similarity, according to the motifs database."""

    # Get the TFs linked by gene ID:
    cmdStr = "select genes.Name from motifs.genes inner join motifs.tfMotifs on motifs.genes.geneID = motifs.tfMotifs.geneID where motifs.tfMotifs.motifPrefix = \"" + str(motifPrefix) + "\";"
    motifsCursor.execute(cmdStr)
    rows = motifsCursor.fetchall()

    tfsByGeneID = map(lambda row: row[0], rows)

    # Get the TFs linked by DBD sequence similarity:
    cmdStr = "select distinct g2.Name from motifs.genes inner join motifs.tfMotifs on motifs.genes.geneID = motifs.tfMotifs.geneID inner join geneProteins on genes.GeneID = geneProteins.GeneID inner join proteins on geneProteins.proteinID = proteins.proteinID inner join dbdBlocks on proteins.proteinID = dbdBlocks.proteinID inner join aaSeq on dbdBlocks.aaSeqID = aaSeq.seqID inner join aaSeqSim on aaSeq.seqID = aaSeqSim.seq1id inner join aaSeq as2 on aaSeqSim.seq2id = as2.seqID inner join dbdBlocks db2 on as2.seqID = db2.aaSeqID inner join proteins p2 on db2.proteinID = p2.proteinID inner join geneProteins gp2 on p2.proteinID = gp2.proteinID inner join genes g2 on gp2.geneID = g2.geneID where motifs.tfMotifs.motifPrefix = \"" + str(motifPrefix) + "\" and aaSeqSim.similarity >= 1;"
    motifsCursor.execute(cmdStr)
    rows = motifsCursor.fetchall()

    tfsBySeqSim = map(lambda row: row[0], rows)

    if (len(rows) > 0):
        tfName = rows[0][0]
    else:
        tfName = None # There are no chip-seq TFs for that motif

    # Fixed bug on July 29th 2013: Uniquify the resulting list before returning:
    tfsForMotifDict = {}
    for tf in tfsByGeneID + tfsBySeqSim:
        tfsForMotifDict[tf] = 1
    return tfsForMotifDict.keys()


def cons2countMatrix(consString, count):
    """Takes a consensus sequence and a count value as input, and generates a
    numpy matrix representing the equivelant pwm. Returns that matrix object
    (rather than a pwm object)."""

    # Set up the alphabet mapping dictionary, in order to map from letters to
    # indeces in the matrix:
    mappingDict = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3}
    aLen = 4

    # Generate a new zeros matrix as long as the alphabet and as wide as the
    # consensus string, to build the count matrix on:
    countMat = numpy.zeros((len(consString), aLen))

    # Set the count value in each column of the matrix, based on the consensus
    # sequence letter at that position, and the specified input count value:
    for seqIdx in range(len(consString)):
        lettIdx = mappingDict[consString[seqIdx]]
        countMat[seqIdx][lettIdx] = count
    return countMat


class pwm(object):
    """A single pwm motif model."""

    def __init__(self, motifData, dataType="memeFile"):
        self.name = None
        self.matrix = None
        self.memeFilePath = None
        self.eValue = None
        self.consensusCount = None # Optional, for count data such as selex
        if (dataType == "memeFile"):
            assert (isinstance(motifData, file))
            # Input data is a MEME file => Initialise accordingly:
            self.initFromMEME(motifData)
            if (self.matrix == []):
                # Matrix was still empty after trying to initialise from
                # MEME file -> Report this as an exception:
                raise ValueError("MEME file had no more motifs in it:" +
                                 motifData.name)
        if (dataType == "xxMotifFile"):
            assert (isinstance(motifData, file))
            # Input data is a xxMotif file => Initialise accordingly:
            self.initFrom_xxMotif(motifData)
            if (self.matrix == []):
                # Matrix was still empty after trying to initialise from
                # xxMotif file -> Report this as an exception:
                raise ValueError("xxMotif file had no more motifs in it:" +
                                 motifData.name)
        elif (dataType == "countsFile"):
            assert (isinstance(motifData, file))
            self.initFromCounts(motifData)
        elif (dataType == "freqMatrix"):
            self.initFromFreqMatrix(motifData)
        elif (dataType == "seqAln"):
            # Input data is a list of strings representing aligned DNA sequences
            # from which this pwm should be constructed:
            self.initFromAlign(motifData)
        elif (dataType == "iniMotifFile"):
            self.initFromInimotifFile(motifData)
        self.logoFilename = None

    def getScoreProf(self, dnaSeq, bgFreqs = [0.25,0.25,0.25,0.25,0.25]):
        """Introduced April 15th 2011.

        Calculates the motif LLR score contribution at each position in the
        specified DNA sequence. Returns the resulting LLR score contributions
        as an array. Sequence positions with the letter "N" result in a zero
        value contribution.

        NOTE: This only works/makes sense when a zero-order background
        model is used. Thus, the bgFreqs is assumed to be an array describing
        such a background model (for A, C, G, T and N)."""

        assert len(dnaSeq) == self.getWidth()

        # Generate a data structure mapping from letter to pwm (and background
        # model) index:
        lett2colIdx = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4, \
                           'a':0, 'c':1, 'g':2, 't':3, 'n':4}

        # The array showing contributions of each letter to the total LLR:
        llrScoreContribs = []

        # Consider each position in the motif...
        for columnIdx in range(self.getWidth()):
            # Current column; bgFreqs[-1] gives the "N" frequency:
            motifColumn = self.getMatrix()[columnIdx] + [bgFreqs[-1]]

            # Get motif likelihood, bg likelihood, and then calculate llr:
            letter = dnaSeq[columnIdx]
            letterColIdx = lett2colIdx[letter]
            motifProb = motifColumn[letterColIdx]
            bgProb = bgFreqs[letterColIdx]
            llrScoreContribs.append(math.log(motifProb/bgProb, 10))

        return llrScoreContribs

    def getName(self):
        return self.name

    def setName(self, name):
        self.name = name

    def getIC_arr(self):
        """Returns an array storing the information content of each column of
        this pwm, in order of motif position."""
        matrix = self.matrix
        ICs = []
        for col in matrix:
            ICs.append(seq_utility.calc_IC(col))
        return ICs

    def getMatrix(self):
        return self.matrix

    def getWidth(self):
        return len(self.matrix)

    def getLogoFilename(self):
        return self.logoFilename

    def getMemeFilePath(self):
        return self.memeFilePath

    def getEValue(self):
        return self.eValue

    def trimLowIC(self, icThresh=0.5, copy=False):
        """Trims off low information-content flanking columns from the motif."""
        trimmedMotif = self.getMatrix()

        # Trim the leading low IC columns...
        while ((len(trimmedMotif) > 0) and
               (seq_utility.calc_IC(trimmedMotif[0]) < icThresh)):
            trimmedMotif = trimmedMotif[1:]

        # Trim the tailing low IC columns...
        while ((len(trimmedMotif) > 0) and
               (seq_utility.calc_IC(trimmedMotif[-1]) < icThresh)):
            trimmedMotif = trimmedMotif[:-1]
        if (not copy):
            # Editing the original motif => Set matrix:
            self.matrix = trimmedMotif
            return None
        else:
            # Making a copy matrix:
            copyMotif = pwm(trimmedMotif, dataType="freqMatrix")
            return copyMotif

    def makeSeqLogo(self, outFilePrefix, format="eps"):
        """Generates an eps sequence logo for this motif, writing it out to the
        specified filename."""

        self.logoFilename = outFilePrefix + "." + format

        tmpMemeFilename = \
            utility.makeTempFilename(outFilePrefix + "_tmp_MEME_file_ceqlogo",
                                     fileSuffix = ".meme")
        tmpMemeFile = open(tmpMemeFilename, 'w')
        self.writeToMEME(tmpMemeFile)
        print >> tmpMemeFile, ""
        tmpMemeFile.flush()
        tmpMemeFile.close()

        # FIXME: Currently, I assume the matrix file exists, and that it is
        # in MEME format. Will need to adapt this in the future when
        # dynamically-generated motifs are run instead:
        seqUtils.run_ceqlogo(tmpMemeFilename, outFilePrefix, format=format)
        cmdlineProgs.deleteFiles([tmpMemeFilename])

    def writeToMAST(self, outFile, pseudo=0.01):
        """Writes the motif out the specified filehandle in MAST format."""

        # The motif matrix data must have been set before this method can
        # be called:
        assert (self.matrix != None)

        lengthStr = str(len(self.matrix))

        hdrString = """MEME version 4.5

ALPHABET= ACGT

strands: + -

Background letter frequencies (from dataset with add-one prior applied):
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF """
        extra = str(self.name) + "\nBL   MOTIF " + str(self.name) + " width= " + lengthStr + " seqs=0\n\nlog-odds matrix: alength= 4 w= " + lengthStr + "\n"
        hdrString = hdrString + extra

        outFile.write(hdrString)

        # Currently just assume uniform background model:
        for col in self.matrix:
            # Add pseudo-counts to the elements of the matrix...
            # FIXME: This could be done better. What pseudo-count
            # value is best? Shouldn't this be in a psfm2llr function?
            newCol = map(lambda prob: (prob+pseudo)/(1+(pseudo*4)), col)

            # From mast source code, it seems the log base is 10, although I'm
            # not 100% sure!:
            try:
                llrCol = map(lambda prob: math.log(prob/0.25, 10), newCol)
            except ValueError, e:
                print >> sys.stderr, "Invalid probability column:", newCol
                raise e
            outFile.write(reduce(lambda p1, p2: str(p1) + " " + str(p2),
                                 llrCol, "") + "\n")
        print >> outFile, ""

    def writeToTRANSFAC(self, outFile, nSeqs=1000, pseudo=0.01):
        """Writes the motif out the specified filehandle in TRANSFAC format."""

        # The motif matrix data must have been set before this method can
        # be called:
        assert (self.matrix != None)

        lengthStr = str(len(self.matrix))

        hdrString = "AC " + self.name + """
XX
TY Motif
ID """ + self.name + """
BF undef
P0\tA\tC\tG\tT\n"""

        outFile.write(hdrString)

        # Currently just assume uniform background model:
        colNum = 0
        for col in self.matrix:
            # Add pseudo-counts to the elements of the matrix...
            newCol = map(lambda prob: (prob+pseudo)/(1+(pseudo*4)), col)
            newColCounts = map(lambda freq: int(freq*1000), newCol)

            colNumStr = ("%2d" % colNum).replace(" ", "0")
            outFile.write(colNumStr +
                          reduce(lambda p1, p2: str(p1) + "\t" + str(p2),
                                 newColCounts, "") + "\n")
            colNum += 1

        print >> outFile, "XX\n//"

    def addPseudoCounts(self, pseudo=0.01):
        """Modifies the psfm for this motif, by adding a specified pseudo count
        at each position."""

        newMatrix = []
        for col in self.matrix:
            newCol = map(lambda prob: (prob+pseudo)/(1+(pseudo*4)), col)
            newMatrix.append(newCol)
        self.matrix = newMatrix

    def writeCountMatrix(self, outFile, nSeqs=100, pseudo=0.01):
        """Writes the motif out as a count matrix, using the specified
        number of sequences."""
        assert (self.matrix != None)

        hdrString = ">" + self.getName()
        outFile.write(hdrString + "\n")
        for col in self.matrix:
            probCol = map(lambda prob: (prob+pseudo)/(1+(pseudo*4)), col)
            countCol = map(lambda prob: int(prob*nSeqs), probCol)
            outFile.write(reduce(lambda p1, p2: str(p1) + " " + str(p2),
                                 countCol, "") + "\n")

    def writeFreqMatrix(self, outFile, nSeqs=100, pseudo=0.01):
        """Writes the motif out as a frequency matrix."""
        assert (self.matrix != None)

        hdrString = ">" + self.getName()
        outFile.write(hdrString + "\n")
        for col in self.matrix:
            probCol = map(lambda prob: (prob+pseudo)/(1+(pseudo*4)), col)
            outFile.write(reduce(lambda p1, p2: str(p1) + " " + str(p2),
                                 probCol, "") + "\n")
        
    def writeToMEME(self, outFile, pseudo=0.01):
        """Writes the motif out the specified filehandle in MEME format."""

        # The motif matrix data must have been set before this method can
        # be called:
        assert (self.matrix != None)

        hdrString = """MEME version 4.5

ALPHABET= ACGT

strands:  + -

Background letter frequencies (from dataset with add-one prior applied):
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF """ + str(self.name) + """
BL   MOTIF """ + str(self.name) + """ width= """ + str(len(self.matrix)) + """ seqs=0
letter-probability matrix: alength= 4 w= """ + str(len(self.matrix)) + " nsites= 10 E= 0\n"

        outFile.write(hdrString)

        # Currently just assume uniform background model:
        for col in self.matrix:
            newCol = map(lambda prob: (prob+pseudo)/(1+(pseudo*4)), col)
            outFile.write(reduce(lambda p1, p2: str(p1) + " " + str(p2),
                                 newCol, "") + "\n")
        outFile.flush()

    def initFromFreqMatrix(self, matrix, name="Unknown"):
        """Initialise the matrix by simply assigning to the specified frequency
        matrix."""

        self.matrix = matrix

    def initFromInimotifFile(self, iniMotifFile, name="Unknown"):
        """Initialise the matrix from an iniMotif file."""

        # Read in the count of the consensus sequence...
        currLine = iniMotifFile.readline()
        elems = currLine.split()
        while ((len(elems) < 2) or (elems[1] != "consensus")):
            currLine = iniMotifFile.readline()
            elems = currLine.split()
        self.consensusCount = int(elems[-1])

        # Skip to the start of the frequency matrix...
        while (currLine[:16] != "Frequency Matrix"):
            currLine = iniMotifFile.readline()

        # Read the matrix information:
        currLine = iniMotifFile.readline()
        freqsTranspose = []
        while (currLine != "\n"):
            elems = currLine.split()
            freqs = map(lambda tok: float(tok), elems[1:])
            freqsTranspose.append(freqs)
            currLine = iniMotifFile.readline()

        freqMatrixTranspose = numpy.matrix(freqsTranspose)
        freqMatrix = freqMatrixTranspose.transpose()
        freqMatrixAsList = freqMatrix.tolist()
        self.matrix = freqMatrixAsList

    def initFromCounts(self, countsInfile, name="Unknown", pseudo=0.01):
        """Initialise the matrix from an input text file of counts, where each
        row has four values specifying the A, C, G, T counts for a given column
        of the psfm."""

        self.name = name

        matrix = []
        for line in countsInfile.readlines():
            elems = line.strip().split()
            countsColumn = map(lambda tok: int(tok), elems)
            countsSum = reduce(lambda count1, count2: count1 + count2,
                               countsColumn)
            freqsColumn = \
                map(lambda count: float(count + pseudo)/float(countsSum + pseudo*4),
                              countsColumn)
            matrix.append(freqsColumn)
        self.matrix = matrix

    def initFromAlign(self, seqAlnList, name="Unknown", pseudo=0.01):
        """Initialise the matrix from a list of aligned DNA sequences."""

        lett_to_idx = \
            {'A':0, 'C':1, 'G':2, 'T':3, \
                 'a':0, 'c':1, 'g':2, 't':3}

        # Initialise count matrix:
        seq_len = len(seqAlnList[0])
        count_matrix = []
        seq_idx = 0
        while (seq_idx < seq_len):
            count_matrix.append([0,0,0,0])
            seq_idx = seq_idx + 1

        for sequence in seqAlnList:
            # Add counts to count_matrix for current sequence:
            seq_idx = 0;
            while (seq_idx < seq_len):
                curr_lett = sequence[seq_idx]
                # If 'N' is found, add a count of 0.25 to each count
                # matrix letter:
                if ((curr_lett == 'N') or (curr_lett == 'n')):
                    count_matrix[seq_idx][0] = count_matrix[seq_idx][0] + 0.25
                    count_matrix[seq_idx][1] = count_matrix[seq_idx][1] + 0.25
                    count_matrix[seq_idx][2] = count_matrix[seq_idx][2] + 0.25
                    count_matrix[seq_idx][3] = count_matrix[seq_idx][3] + 0.25
                else:
                    lett_idx = lett_to_idx[curr_lett]
                    count_matrix[seq_idx][lett_idx] = \
                        count_matrix[seq_idx][lett_idx] + 1
                seq_idx = seq_idx + 1
                self.name = name

        freqMatrix = []
        for countsColumn in count_matrix:
            countsSum = reduce(lambda count1, count2: count1 + count2,
                               countsColumn)
            freqsColumn = \
                map(lambda count: float(count + pseudo)/float(countsSum + pseudo*4),
                              countsColumn)
            freqMatrix.append(freqsColumn)
        self.matrix = freqMatrix

    def initFrom_xxMotif(self, xxMotif_infile):
        # Store name of file as attribute of motif:

        self.xxMotifFilePath = xxMotif_infile.name

        # Parse a single motif from the assumed open file...

        # Throw away lines until the next motif line is reached...
        currLine = xxMotif_infile.readline()
        if (currLine == ""):
            raise ValueError(xxMotif_infile.name + " contains no more motifs!")

        while ((currLine != "") and (currLine[:5] != "Motif")):
            currLine = xxMotif_infile.readline()

        if (currLine == ""):
            raise ValueError(xxMotif_infile.name + " contains no more motifs!")

        # Parse the motif line, to obtain the motif name:
        elems = currLine.strip().split()
        motifName = elems[1][:-1]
        motifEvalue = float(elems[-1])
        self.name = motifName
        self.eValue = motifEvalue

        # Parse the frequency matrix:
        freqMatrixList = []
        for lettIdx in [1,2,3,4]:
            currLine = xxMotif_infile.readline()
            elems = currLine.strip().split()
            freqsCol = map(lambda tok: float(tok), elems[1:])
            freqMatrixList.append(freqsCol)

        freqMatrix = numpy.matrix(freqMatrixList)
        freqMatrixTranspose = freqMatrix.transpose()
        freqMatrixTransposeList = freqMatrixTranspose.tolist()
        self.matrix = freqMatrixTransposeList

    def initFromMEME(self, meme_infile):
        # Initialise this motif from the specified motifData input file...

        # Store the name of the file as an attribute of this motif:
        self.memeFilePath = meme_infile.name

        # Parse a single motif from an open MEME file...

        # While the currline doesn't start with "MOTIF", discard the line and
        # get the next...
        currLine = meme_infile.readline()
        if (currLine == ""):
            raise ValueError(meme_infile.name + " contains no more motifs!")

        # Ignore lines until "MOTIF" is reached:
        while ((currLine != "") and (currLine[:5] != "MOTIF")):
            currLine = meme_infile.readline()

        if (currLine == ""):
            raise ValueError(meme_infile.name + " contains no more motifs!")

        # Parse the motif line, to obtain the motif name:
        elems = currLine.strip().split()
        motifName = elems[1]

        # Discard lines until the  "letter-probability" line is reached...
        while ((currLine != "") and (currLine[:18] != "letter-probability")):
            currLine = meme_infile.readline()

        if (currLine == ""):
            raise ValueError(meme_infile.name + " contains no matrix info" + \
                                 " for matrix " + motifName)

        # Parse the motif stats; extract and save the E-value of the motif:
        elems = currLine.strip().split()
        eValue = float(elems[9])
        nsites = int(elems[7])

        # Instantiate the motif with the parsed values:
        self.name = motifName
        self.eValue = eValue
        self.nsites = nsites

        # Parse the actual matrix data...
        currLine = meme_infile.readline()
        matrix = []
        # Modified April 17th 2013; allowing new line character to indicate
        # end of current motif:
        while ((currLine != "") and (currLine[:3] != "---") and (currLine != "\n")):
            elems = currLine.strip().split()
            matrixCol = map(lambda tok: float(tok), elems)
            matrix.append(matrixCol)
            currLine = meme_infile.readline()
        self.matrix = matrix

    def getRevCompMatrix(self):
        """Returns a raw frequency matrix representing the reverse complement of this matrix."""
        matrixCopy = copy.deepcopy(self.matrix)
        matrixCopy.reverse()
        for col in matrixCopy:
            col.reverse()
        return matrixCopy

    def getConsensus(self):
        """Returns the consensus sequence of this motif."""

        letters = ["A","C","G","T"]

        maxPositions = numpy.array(self.matrix).argmax(1)
        consensus = reduce(lambda lett1, lett2: lett1+lett2,
                           map(lambda maxPos:letters[maxPos], maxPositions))
        return consensus

    def loadIntoDb(self, analysisID, cycle, db):
        """Loads this motif into the specified SELEX database, with the
        specified analysisID and SELEX cycle number that it was generated
        from."""

        # Obtain the ensemble ID of this factor, using the analysisID and a
        # mysql query...
        
        # Find the Barcode and Batch studied with the specified analysis:
        cursor=db.cursor()
        cmdStr = "select BarcodeAnalysed, BarcodeBatchAnalysed from SELEX_Analyses where (AnalysisID = " + str(analysisID) + ");"
        cursor.execute(cmdStr)
        rows = cursor.fetchall()

        # Return the first element of the first result from the query run:
        (barcode, batch) = rows[0]

        # Obtain the ensembleID, by using a SelexSample object:
        selexSample = selexDb.SelexSample(barcode, batch, cycle)
        ensembleID = selexSample.getEnsembleID()

        # Obtain the width of the motif:
        motWidth=self.getWidth()

        # Obtain the consensus sequence:
        consensus=self.getConsensus()
        consensusCount=self.consensusCount

        # Insert a new entry into the Motifs table using a new mysql command,
        # inserting analysisID, ensembleID, width, cycle, consensus sequence
        # and consensus count:
        if (ensembleID != None):
            cmdStr = "INSERT INTO Motifs (AnalysisID, EnsembleID, Width, SelexCycle, ConsensusSequence, ConsensusCount) VALUES (" + str(analysisID) + ",\"" + str(ensembleID) + "\"," + str(motWidth) + "," + str(cycle) + ",\"" + consensus  + "\"," + str(consensusCount) + ");"
        else:
            cmdStr = "INSERT INTO Motifs (AnalysisID, Width, SelexCycle, ConsensusSequence, ConsensusCount) VALUES (" + str(analysisID) + "," + str(motWidth) + "," + str(cycle) + ",\"" + consensus  + "\"," + str(consensusCount) + ");"
    
        # Run that command:
        cursor=db.cursor()
        cursor.execute(cmdStr)

        # Obtain the new autoincremented motifID value:
        newMotifID = db.insert_id()

        # Insert data for the columns...
        for columnIdx in range(len(self.matrix)):
            currColumn = self.matrix[columnIdx]

            # Calculate information content of that column:
            ic = seq_utility.calc_IC(currColumn)

            # Insert a new entry into MotifColumns, with the above motifID,
            # columnIdx, information content, and A, C, G and T frequencies:
            cmdStr = "INSERT INTO MotifColumns (MotifID, ColumnNumber, InformationContent, A, C, G, T) VALUES (" + str(newMotifID) + "," + str(columnIdx) + "," + str(ic) + "," + str(currColumn[0]) + "," + str(currColumn[1]) + "," + str(currColumn[2]) + "," + str(currColumn[3]) + ");"
    
            # Run that command:
            cursor=db.cursor()
            cursor.execute(cmdStr)

    def permuteMatrix(self):
        """Permutes the columns of this matrix."""
        random.shuffle(self.matrix)

    def isCrap(self):
        """Returns True if this motif is dodgy, False otherwise."""

        # A bit tricky; Get a trimmed version of the motif, and then compute
        # whether the motif is crap based on the trimmed version and the
        # original version...

        trimmedMot = self.trimLowIC(icThresh=0.5, copy=True)
        if (self.eValue > 0.01):
            return True
        if (trimmedMot.getWidth() <= 10):
            return False
        if (self.eValue < 10**-30):
            return False
        motICs = trimmedMot.getIC_arr()
        letters = trimmedMot.getConsensus()
        lettersAtHighIC = {}
        columnIdx = 0
        while (columnIdx < trimmedMot.getWidth()):
            currLetter = letters[columnIdx]
            ic = motICs[columnIdx]
            if (ic >= 1.5):
                if (not lettersAtHighIC.has_key(currLetter)):
                    lettersAtHighIC[currLetter] = 1
            columnIdx += 1
        if (len(lettersAtHighIC.keys()) >= 3):
            return False
        return True


class JASPAR_pwm(pwm):
    """A pwm derived from a high-throughput SELEX experiment."""
    def __init__(self, db, motifID):
        super(JASPAR_pwm, self).__init__(None, dataType="")

        # Obtain the JASPAR ID and collection values for this motifID:
        cmdStr = "SELECT COLLECTION, BASE_ID, NAME, VERSION from MATRIX where ID = " + str(motifID) + ";"
        cursor=db.cursor()
        cursor.execute(cmdStr)
        rows = cursor.fetchall()

        assert(len(rows) > 0) # Otherwise motifID is an invalid input.

        self.collection = rows[0][0]
        self.jasparID = rows[0][1]
        tfName = rows[0][2]
        version = rows[0][3]

        # Set the name of the matrix to the jaspar code and tf name:
        self.setName(self.jasparID + "_" + str(version) + "_" + tfName)

        # Obtain the motif data:
        cmdStr = "SELECT col, row, val from MATRIX_DATA where ID = \"" + str(motifID) + "\";"
        cursor=db.cursor()
        cursor.execute(cmdStr)
        rows = cursor.fetchall()

        positions = map(lambda tup: tup[0], rows)
        maxColNum = max(positions)

        # Parse the output into the matrix attribute...

        # Obtain the count matrix:
        countMatrix = numpy.zeros((maxColNum, 4))
        lett2col = {'A':0, 'C':1, 'G':2, 'T':3}
        for tup in rows:
            row = tup[0] - 1
            col = lett2col[tup[1]]
            countMatrix[row][col] = tup[2]

        # Divide all entries by row sums to obtain freq matrix:
        freqMatrix = numpy.zeros((maxColNum, 4))
        for rowIdx in range(maxColNum):
            rowSum = float(sum(countMatrix[rowIdx]))
            for colIdx in range(4):
                freqMatrix[rowIdx][colIdx] = \
                    (countMatrix[rowIdx][colIdx])/rowSum
        
        self.matrix = freqMatrix.tolist()


class SELEX_pwm(pwm):
    """A pwm derived from a high-throughput SELEX experiment."""
    def __init__(self, db, motifID):
        super(SELEX_pwm, self).__init__(None, dataType="")

        # Obtain the cycle and consensusCount values for this motifID:
        cmdStr = "SELECT SelexCycle, ConsensusCount from Motifs where MotifID = \"" + str(motifID) + "\";"
        cursor=db.cursor()
        cursor.execute(cmdStr)
        rows = cursor.fetchall()

        self.cycle = rows[0][0]
        self.consensusCount = rows[0][1]

        # Obtain the motif data:
        cmdStr = "SELECT ColumnNumber, A, C, G, T from MotifColumns where MotifID = \"" + str(motifID) + "\";"
        cursor=db.cursor()
        cursor.execute(cmdStr)
        rows = cursor.fetchall()
        matrixData = list(rows)
        matrixData.sort(cmpTups1)

        # Obtain a name for the motif:
        cmdStr = "select tf.HGNC_Name, \"w\", m.Width, \"cyc\", m.SelexCycle, sa.BarcodeBatchAnalysed, sa.BarcodeAnalysed from Motifs m inner join SELEX_Analyses sa on m.AnalysisID = sa.AnalysisID inner join TranscriptionFactors tf on m.EnsembleID = tf.EnsembleID where m.MotifID = " + str(motifID) + ";"
        cursor=db.cursor()
        cursor.execute(cmdStr)

        rows = cursor.fetchall()
        name = reduce(lambda tok1, tok2: str(tok1) + "_" + str(tok2), rows[0])

        self.setName(name)

        # Parse the output into the matrix attribute:
        matrix = map(lambda col: col[1:], matrixData)
        self.matrix = matrix

    def getConsCount(self):
        return self.consensusCount


def cmpTups1(tup1, tup2):
    if (tup1[0] < tup2[0]):
        return -1
    elif (tup1[0] == tup2[0]):
        return 0
    else:
        return 1


class motifLibrary:
    """A library of motifs. Introduced on 9th November 2010 in order to allow
    MEME-derived motifs to be added to an existing set of motifs more
    seamlessly."""

    def __init__(self, inputLibrary):
        self.motifs = {}
        if (isinstance(inputLibrary, list)):
            # Library specifies name of file containing motif locations =>
            # initialise as such:
            self.initFromMEME(inputLibrary)

    def getMotif(self, motifName):
        return self.motifs[motifName]

    def getMotifs(self):
        return self.motifs

    def initFromMEME(self, motifLibraryList):
        """motifLibraryList is a list of strings specifying the absolute
        paths of motifs in the input library."""

        # Set up the motif library from the input file...

        # For each motif filename in the input library...
        for currMotifFilename in motifLibraryList:
            if (not os.path.exists(currMotifFilename)):
                # Don't quit if a motif file isn't found; just skip it:
                print >> sys.stderr, "WARNING: Specified motif file \"" + \
                    currMotifFilename + "\" does not exist."
            else:
                # Attempt to generate a new motif object from that file, using
                # the motif constructor:
                try:
                    currMotif = pwm(open(currMotifFilename))
                    currMotifName = currMotif.getName()
                    if (self.motifs.has_key(currMotifName)):
                        print >> sys.stderr, "WARNING: motif \"" + \
                            currMotif.getName() + \
                            "\" was encountered more than once." + \
                            " Excluding second instance from the input library."
                    else:
                        # Non-redundant motif was successfully parsed => Add
                        # it to the library of motifs:
                        self.motifs[currMotifName] = currMotif
                except ValueError, e:
                    print >> sys.stderr, "WARNING: ValueError occured \"" + \
                        "whilst parsing motif file \"" + currMotifFilename + \
                        "\". Excluding that motif from the input library."
                    print >> sys.stderr, e

    def addMotif(self, motif):
        """Add a single specified motif to this library."""
        assert(not self.motifs.has_key(motif.getName()))
        self.motifs[motif.getName()] = motif


class scannedReg:
    """This class encapsulates the general idea of a genomic region that has
    been scanned by a set of motifs."""

    def __init__(self, region):
        # Store the attributes (in particular the genomic coordinates) of the
        # genomic region:
        self.reg = region # A BED_line object.
        
        # Set up an intially-empty dictionary of motif hits. It will be
        # populated with motif name as keys and motifOccurence as values:
        self.hits = {}

    def addHit(self, motifHit):
        """Stores the specified motifOccurence object under it's name.
        Throws an Exception if a motif by that name has already been
        registered for this sequence."""

        if (self.hits.has_key(motifHit.getMotif().getName())):
            raise Exception("Region already had motif instance of that name.")
        else:
            self.hits[motifHit.getMotif().getName()] = motifHit

    def getHit(self, motifName):
        """Returns the occurence of the specified motif for this region.
        Returns None if no hits are registered for the region."""
        if (self.hits.has_key(motifName)):
            return self.hits[motifName]
        else:
            return None

    def getDisp(self, mot1_name, mot2_name):
        """Returns a (displacement, sameStrand) tuple showing the
        displacement from motif1 to motif2 and whether they occur on the
        same strand of this sequence. Returns None if the two motifs overlap.
        Throws a KeyError exception if either motif is not present in the
        sequence."""
        mot1_hit = self.hits[mot1_name]
        mot2_hit = self.hits[mot2_name]
        disp1to2 = mot1_hit.get_disp(mot2_hit)

    def getBED(self):
        """Returns the BED_line object showing the coordinates of this
        region."""
        return self.reg

    def getMotCoords(self, motifName, flankWidth=0):
        """Returns a strand-specific BED_line object representing the
        *genomic* coordinates of the occurrence of the specific motif in this
        scanned region. Returns None if no motif hit exists in this region for
        the specified motif."""

        motifHit = self.getHit(motifName)
        if (motifHit == None):
            return None

        genomicHitLoc = chipseq_analysis.bestHit_scanner.getAbsHitCoords(self.getBED(), motifHit, flankWidth=flankWidth)
        return genomicHitLoc
