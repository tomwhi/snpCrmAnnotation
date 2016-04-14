#!/usr/bin/env python

import commands, cmdlineProgs, fileinput, os, os.path, pdb, random, re, sys, string, time
#import pygraph.algorithms.accessibility as pygraphAc
import matplotlib.pyplot as plt
import matplotlib
import numpy, scipy.stats


def listFtpDir(ftpDir):
    """Retrieves the contents of a specified ftp directory.
    Returns an array of the content filenames"""

    # Current uses a command-line text web browser to accomplish this,
    # since I couldn't find a way to use ftp without having to start
    # an interactive shell.

    # Sleep for 1 second, as a precaution against bombarding the ftp
    # server with requests:
    time.sleep(1)

    # Retrieve the lynx output:
    tmpFtpContents = makeTempFilename("tmpFtpContents")
    cmdStr = "lynx -dump " + ftpDir + " > " + tmpFtpContents
    cmdlineProgs.runCommand(cmdStr, True)

    # Extract all filenames:
    allText = "\n".join(open(tmpFtpContents).readlines())
    listedFiles = re.findall("ftp://[^\n]+", allText)

    cmdlineProgs.deleteFiles([tmpFtpContents])

    # Return the filenames:
    return listedFiles


def fileSizeComp(file1Loc, file2Loc):
    """Comparator function that compares two file locations based on their
    disk usage."""

    file1diskUsage = int(commands.getstatusoutput("du " + file1Loc)[1].split("\t")[0])
    file2diskUsage = int(commands.getstatusoutput("du " + file2Loc)[1].split("\t")[0])
    return file1diskUsage - file2diskUsage


def getFileDir(fileLoc):
    """Returns path to directory of specified file."""
    lastBackslashLoc = fileLoc.rfind("/")
    chipAlnDir = fileLoc[:lastBackslashLoc]

    return chipAlnDir


def is_int(token):
    """Returns true if the string represents and integer, false otherwise."""
    try:
        int(token)
        return True
    except ValueError:
        return False

def pgv2adjMatrix(pgvGraph, outFile):
    """Writes the specified pgv graph out as an adjacency matrix."""
    # Print nodes header and adjtable header:
    nodesLine = "set NODES := "
    for node in pgvGraph.nodes():
        nodesLine = nodesLine + " \"" + str(node) + "\""
    print >> outFile, nodesLine + ";"
    print >> outFile, "param AdjTable:" + reduce(lambda tok1, tok2: tok1 + " " + tok2, map(lambda nodeName: "\"" + nodeName + "\"", pgvGraph.nodes())) + ":="

    # Print the adjacency matrix body:
    for node1 in pgvGraph.nodes():
        adjRow = ""
        for node2 in pgvGraph.nodes():
            if pgvGraph.has_edge(node1, node2):
                adjRow = adjRow + " 1"
            else:
                adjRow = adjRow + " 0"
        outLine = "\"" + node1 + "\"" + adjRow
        print >> outFile, outLine
    print >> outFile, ";"

    
def findConnComps(graph):
    """Finds all connected components in the specified graphviz graph object,
    and stores them as subgraphs in the original graph."""

    # Will use the module pygraph to facilitate this...

    # Use pygraph to find a component label for each node in the graph:
    try:
        node2component = \
            pygraphAc.connected_components(graph)

        # Generate a set of components, represented as a dictionary where
        # each key is a component label, and each value is a node...
        component2nodes = {}
        for node in graph.nodes():
            componentLabel = node2component[str(node)]
            if (component2nodes.has_key(componentLabel)):
                # There is an existing component recognised already => Record
                # the node in that one:
                component2nodes[componentLabel].append(node)
            else:
                # Start a new component list containing only this node:
                component2nodes[componentLabel] = [node]

        # Generate a subgraph for each of those lists, and register it in
        # the input graph using it's "subgraph()" method...
        for nodeList in component2nodes.values():
            graph.subgraph(nodeList)

    except KeyError, e:
        # Pygraphviz package seems to have some kind of horrible intermittent
        # bug. To stop this from crashing my program altogether, I will
        # return without identifying connected components if this occurs.
        print >> sys.stderr, "WARNING: pygraphviz bug seems to have occurred."+\
            " Returning without finding connected components."
        return


def getRelPath(startPath, endPath):
    """Takes a "starting point" directory and an "end point" file (or directory) as input. Calculates the relative path from the specified starting point to the specified end point location. Returns the normalised path."""

    # Calculate how many steps "up" are need to get from the specified
    # startPath to the *shared* section of the specified endPath...
    startPathElems = startPath.strip("/").split("/")
    endPathElems = endPath.strip("/").split("/")

    nSharedElems = 0
    minLen = min(len(startPathElems), len(endPathElems))
    elemIdx = 0
    if (minLen > 0):
        currStartPathElem = startPathElems[elemIdx]
        currEndPathElem = endPathElems[elemIdx]
        while ((elemIdx < minLen) and (currStartPathElem == currEndPathElem)):
            # The current directory path elements match => record this...
            nSharedElems = nSharedElems + 1
            elemIdx = elemIdx + 1

            if (elemIdx < minLen):
                currStartPathElem = startPathElems[elemIdx]
                currEndPathElem = endPathElems[elemIdx]

    # Number of steps "up" required is the number of steps need to get to
    # root, minus the number of shared elements:
    startToShared = "../"* (len(startPath.strip("/").split("/")) - nSharedElems)

    # Reconstitute the remainder of the endPath, from the shared section
    # onwards...
    sharedToEnd = reduce(lambda elem1, elem2: elem1 + "/" + elem2,
                         endPathElems[nSharedElems:])

    relPath = startToShared + sharedToEnd
    return os.path.normpath(relPath)


def addHtmlBreaks(inStr, nChars):
    """Adds an html line-break every nth character, and returns the resulting
    string."""
    strWithBreaks = ""
    strIdx = 0
    while (strIdx < len(inStr)):
        if (((strIdx + 1) % nChars) == 0):
            strWithBreaks = strWithBreaks + "</br>"
        strWithBreaks = strWithBreaks + inStr[strIdx]
        strIdx = strIdx + 1
    return strWithBreaks


def listSum(x, y):
    if (len(x) != len(y)):
        print >> sys.stderr, "ERROR: mismatch in lengths"
        sys.exit(1)
    sumList = []
    idx = 0
    while (idx < len(x)):
        sumList.append(x[idx] + y[idx])
        idx = idx + 1
    return sumList

def getCommandUsed():
    """Returns the command-line contents used."""
    return reduce(lambda str1, str2: str1 + " " + str2, sys.argv)


class FormatError(Exception):
    def __init__(self, msg):
        self.args = msg
    def __str___(self):
        return repr(self.args)


def makeTempFilename(filePrefix, fileSuffix = ".tmp", check=True, tmpDir="/tmp"):
    """Generates a new temporary filename with the specified prefix, and the
    current process id in the filename. Raises an exception if the specified
    filename already exists. "check" specifies whether to check if the
    specified file already exists.

    August 5th 2013: Changed the default behaviour. If "check" is True, then
    now it will by default generate a different unique name, while the currently-generated
    name already exists (to a maximum of 10 attempts).
"""

    # Generate a temporary filename:
    pid = os.getpgid(0)
    filenameEnd = fileSuffix
    #filenameStartLen = 254 - len(filenameEnd)
    #filenameStart = filePrefix[:filenameStartLen] + "_" + str(pid)
    filenameStart = filePrefix + "_" + str(pid)

    ranAttemptCount = 0
    tmp_filename = tmpDir + "/" + filenameStart + "_" + str(ranAttemptCount) + "_" + fileSuffix

    # Modified in August 5th 2013:
    if check:
        # Checking if the file exists. Now (new default behaviour), try to
        # generate a unique name, up to 10 times before giving up...
        while (os.access(tmp_filename, os.F_OK) and ranAttemptCount < 10):
            ranAttemptCount += 1
            tmp_filename = filenameStart + "_" + str(ranAttemptCount) + "_" + fileSuffix

        if (os.access(tmp_filename, os.F_OK)):
            # Even after trying several times to generate a unique temporary
            # filename, the processed temporary filename is invalid, since a file
            # of that name already exists => Abort with an error:
            print >> sys.stderr, "ERROR: Proposed temporary file, " + \
                tmp_filename + " already exists."
            raise Exception()

    return tmp_filename


def makeDir(dirName):
    """Makes a new directory with the specified name. If the directory already
    exists then raise a new exception."""

    if (os.access(dirName, os.F_OK)):
        raise Exception("Directory already exists: " + dirName)

    return commands.getstatusoutput("mkdir " + dirName)


def copyFiles(fileNames, targetDir):
    """Calls linux cp."""

    filesStr = reduce(lambda name1, name2: name1 + " " + name2, fileNames)
    cmdStr = "cp " + filesStr + " " + targetDir
    return commands.getstatusoutput(cmdStr)


def ksTest():
    return # XXX Implement later, or perhaps get rpy working instead.


def binomial(s, n, p):
    """Returns the probability of at least <s> successes in <n> 
	independent Bernouli trails each with probability
	of success <p>."""

    return scipy.stats.binom.sf(s-1, n, p)

    # cmdStr = "binomial " + str(s) + " " + str(n) + " " + str(p)
    # try:
    #     #print >> sys.stderr, "TRACE1:", scipy.stats.binom.sf(s-1, n, p)
    #     #print >> sys.stderr, "TRACE2:", float(commands.getoutput(cmdStr).split()[-1])
    #     return float(commands.getoutput(cmdStr).split()[-1])
    # except ValueError, e:
    #     print >> sys.stderr, cmdStr
    #     raise(e)


class R_plotter:
    """A class for generating plots using R. I won't implement too much here
    because I might be reinventing the wheel (I believe there is an "Rpy"
    package)."""

    def __init__(self):
        pass

    def makeHist(self, outputFilePrefix, vals):
        """Generates a graphics file showing a histogram of the input values."""

        # Currently I will implement this in R. A disadvantage of this is
        # that it must be run in an X-term on MAC OS X.

        tmpValsFilename = makeTempFilename("ValsFile")
        tmpValsFile = open(tmpValsFilename, 'w')
        for val in vals:
            if (val != None):
                print >> tmpValsFile, val
                tmpValsFile.flush()
                tmpValsFile.close()

        histCmd = """png('""" + outputFilePrefix + """.png')
        hist(scan('""" + tmpValsFilename + """'))
        dev.off()
        """

        tmp_R_filename = makeTempFilename("R_cmds")
        tmp_R_file = open(tmp_R_filename, 'w')
        print >> tmp_R_file, histCmd
        tmp_R_file.flush()
        tmp_R_file.close()

        # Run R on the specified file, thus (hopefully) generating the required
        # histogram:
        cmdresult = commands.getstatusoutput("R --no-save < " + tmp_R_filename)
        print >> sys.stderr, "Result of R histogram commands:", cmdresult
        
        deleteFiles([tmp_R_filename, tmpValsFilename])


def hammingDist(str1, str2):
    """Returns the hamming distance between the two equal-length strings."""
    assert (len(str1) == len(str2))
    hDist = 0
    for idx in range(len(str1)):
        char1 = str1[idx]
        char2 = str2[idx]
        if (char1 != char2):
            hDist = hDist + 1
    return hDist


def strDist(shortStr, longStr, maxOverhang):
    """Calculates the distance between two strings, allowing the specified shift
    that states how much the shorter string is allowed to overhang the longer
    string when generating the alignment between them."""

    # Figure out the possible displacement values of the start of the short
    # string relative to the start of the long string, allowing the specified
    # maximum amount of overhang:
    minDisp = -maxOverhang
    maxDisp = len(longStr) - len(shortStr) + maxOverhang

    # Minimum distance observed thus far; infinity:
    minDist = len(longStr) + 1000000

    # Consider each possible alignment of the two strings...
    for disp in range(minDisp, maxDisp + 1):
        # Grab the two sequence sections that overlap in the current alignment:
        if (disp < 0):
            # Start of short string overhangs start of longer string:
            str1 = shortStr[-disp:]
            str2 = longStr[:len(str1)]
        else:
            # Start of short string is past the start of the longer string:
            str1 = shortStr[:len(longStr) - disp]
            str2 = longStr[disp:(disp+len(str1))]

        # Calculate distance value of the current alignment. Each overhanging
        # letter (on either side) counts the same as a single mismatch...
        nOverHangs = abs(disp) + abs((len(shortStr) + disp) - len(longStr))
        currDist = abs(nOverHangs) + hammingDist(str1, str2)
#        pdb.set_trace()
        if (currDist < minDist):
            minDist = currDist
    return minDist


def cmpTup2(tupA, tupB):
    """A comparator function that compares two tuples on the basis of the
    value of their second element."""
    if (tupA[1] < tupB[1]):
        return -1
    elif (tupA[1] == tupB[1]):
        return 0
    else:
        return 1


def calc_pi_1(pValsFilename):
    """Runs a q-value analysis on an input file of p-values, and reports the
    resulting pi_1 value."""

    # FIXME: Need to introduce some way of checking whether the command is
    # successful, and reporting an error otherwise.

    outputFilename = pValsFilename + ".pi0" # Hack
    cmd = "qvalue --pi-zero " + pValsFilename + " > " + \
        outputFilename

    print >> sys.stderr, "Running qValue, with command:\n", cmd
    cmdResult = commands.getstatusoutput(cmd)
    print >> sys.stderr, "Command result:", cmdResult

    try:
        outStr = cmdResult[1].split("\n")[1]
        piZero = float(outStr.strip().split("=")[1])
    except Exception, e:
        print >> sys.stderr, "WARNING: Could not extract pi-zero from file " + \
            outputFilename
        print >> sys.stderr, e
        return -1
    return 1 - piZero


def addPlotToAxis(matrix, plotAxes, colour="black"):
    """Adds a line plot to the specified matplotlib "axis" object,
    for the specified 2D matrix of data."""

    npMatrix = numpy.matrix(matrix)
    matTranspose = npMatrix.transpose()
    xVals = matTranspose[0].tolist()[0]
    yVals = matTranspose[1].tolist()[0]

    # Generate plot without error bars:
    barPlot = plotAxes.plot(xVals, yVals, linestyle="dashed", marker="o",
                            color=colour)

    if len(matTranspose) > 2:
        # Matrix includes error bar info:
        assert len(matTranspose) == 4

        errBarsLower = matTranspose[2].tolist()[0]
        errBarsUpper = matTranspose[3].tolist()[0]
        errBarsCombined = [errBarsLower, errBarsUpper]

        # Determine how often to plot the error bars; we want 20 in total:
        nBars = 20
        errBarSpacing = len(errBarsLower)/nBars

        # Generate error bars:
        plotAxes.errorbar(xVals, yVals, color=colour, errorevery=errBarSpacing,
                          yerr=errBarsCombined)


def makePlot(matrix, outputPlotFile, xlabel="No x-axis label specified",
             ylabel="No y-axis label specified",
             title = "No title specified", plotFormat="png"):
    """Generates a (scatter)plot for the specified matrix of data, to the specified
    output file. If the matrix has two columns, then generate without
    error-bars. Otherwise, it must have four columns, and the last two
    columns specify the upper and lower bounds of the error bars,
    respectively."""

    plt.cla()
    plotAxes = plt.axes()
    addPlotToAxis(matrix, plotAxes)
    plotAxes.xaxis.set_label_text(xlabel)
    plotAxes.yaxis.set_label_text(ylabel)
    plotAxes.set_title(title)
    plt.savefig(outputPlotFile, format=plotFormat)


def makeHeatmap(inputData, outFilePrefix, yLabelOrder=None, xLabels=None,
                cmap=plt.cm.OrRd):
    """Generate a heatmap from the input dictionary (specifying Y-axis labels
    and rows of Z-axis values), and save the figure to the specified file."""

    if (yLabelOrder == None):
        yTickLabels = inputData.keys()
        yTickLabels.sort()
    else:
        assert(len(yLabelOrder) == len(inputData.keys()))
        yTickLabels = yLabelOrder

    # NOTE: Have to reverse the order as pcolor counts from the bottom up:
    yTickLabels.reverse()

    # Generate a numpy 2D array from the input dictionary...

    # Generate an empty matrixAsList:
    matrixAsList = []

    # For each row of the heatmap specified as an entry in the dictionary,
    # in the order specified by the optional array of Y-axis labels (or just
    # the dictionary order, if the optional array is None)...
    for yLabel in yTickLabels:
        # Append the row to the matrixAsList:
        matrixAsList.append(inputData[yLabel])

    # Convert the resulting matrixAsList to a numpy array:
    matrixAsArray = numpy.array(matrixAsList)

    # Generate bounding boxes to be input into pcolor...

    # Generate an xBounds array...
    xBounds = [0]
    xMatrixIdx = 1
    while (xMatrixIdx <= len(inputData.values()[0])):
        xBounds.append(xMatrixIdx)
        xMatrixIdx = xMatrixIdx + 1

    # Generate a yBounds array...
    yBounds = [0]
    yMatrixIdx = 1
    while (yMatrixIdx <= len(inputData.keys())):
        yBounds.append(yMatrixIdx)
        yMatrixIdx = yMatrixIdx + 1

    # Generate the plot, using matplotlib.pcolor...
    plt.cla()
    
    plt.subplots_adjust(left=0.3)
    print >> sys.stderr, "Running pylab.pcolor command..."
    plt.pcolor(numpy.array(xBounds), numpy.array(yBounds), matrixAsArray,
             cmap=cmap)
#    plt.colorbar()
    print >> sys.stderr, "Finished running command."
    # Calculate fontsize of ytick labels dynamically, according to the
    # number of rows being plotted:
    fontsize=4 # FIXME: Implement this here...
    ax = plt.gca()
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    # Make the border on the left large and the labels small:
#    plt.axes([.2,.15,.95-0.2,.95-0.2])
#    plt.rcParams.update({'ytick.labelsize':"small"})

    # Generate axis label ticks...

    # Generate Y labels...
    yTickLocs = []
    currYlabelLoc = 0.5
    for yAxisLabel in yTickLabels:
        yTickLocs.append(currYlabelLoc)
        currYlabelLoc = currYlabelLoc + 1
    plt.yticks(numpy.array(yTickLocs), yTickLabels)

    # Generate X labels...

    # Will just mark the first, middle, and last columns...
    firstXlab_idx = 0
    middleXlab_idx = round((len(matrixAsArray[0]) - 1)/2.0)
    lastXlab_idx = len(matrixAsArray[0]) - 1

    firstXlab_loc = firstXlab_idx + 0.5
    middleXlab_loc = middleXlab_idx + 0.5
    lastXlab_loc = lastXlab_idx + 0.5

    xTickLocs = [firstXlab_loc, middleXlab_loc, lastXlab_loc]
    if (xLabels == None):
        xTickLabels = [firstXlab_idx, middleXlab_idx, lastXlab_idx]
    else:
        xTickLabels = xLabels
    plt.xticks(numpy.array(xTickLocs), xTickLabels)

    print >> sys.stderr, "Saving figure..."
    plt.savefig(outFilePrefix + ".png", format="png", bbox_inches='tight', dpi=200)
    print >> sys.stderr, "Saved."


def sampleLines(inputFilename, outFile, fractionSampled, commentChar="#"):
    """Samples lines from the specified text file, outputting to the specified
    file."""

    # August 13th 2013: Moved from sampleLines.py, to facilitate re-use in
    # peak-calling saturation analysis.

    # Find out the number of (non-comment) lines in the file:
    nonCommentLineIdxs = []
    lineIdx = 0
    for line in open(inputFilename).xreadlines():
        if (len(line.strip()) == 0):
            # The line is not a comment => Count it:
            nonCommentLineIdxs.append(lineIdx)
        elif (line[0] != commentChar):
            # The line is not a comment => Count it:
            nonCommentLineIdxs.append(lineIdx)
        lineIdx += 1
    #print >> sys.stderr, "Number of non-comment lines in file:", \
    #    len(nonCommentLineIdxs)

    # Calculate the number of items to randomly sample:
    nSampled = int(len(nonCommentLineIdxs)*fractionSampled)

    #print >> sys.stderr, "Number of lines being sampled:", nSampled

    # Sample that many line numbers; generate an array of such line numbers:
    sampledLineNumbers = random.sample(nonCommentLineIdxs, nSampled)
    sampledLineNumbers.sort()

    # Proceed through the file and print out a line whenever it has one of
    # the sampled line numbers, or if it is a comment line:
    lineIdx = 0
    currSampleIdx = 0
    try:
        for line in open(inputFilename).xreadlines():
            if (currSampleIdx < len(sampledLineNumbers)):
                nextLineIdxToSample = sampledLineNumbers[currSampleIdx]
            else:
                nextLineIdxToSample = -1
            if (lineIdx == nextLineIdxToSample):
                # The current line is to be sampled => print it:
                print >> outFile, line.strip()
                currSampleIdx += 1
            elif (len(line.strip()) > 0 and line[0] == commentChar):
                # The current line is a comment line => Print it:
                print >> outFile, line.strip()
            lineIdx += 1
    except Exception, e:
        pdb.set_trace()
        x = 1

