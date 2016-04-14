#!/usr/bin/env python

import commands, os, sys, utility

"""A module for running various command-line programs in python."""


def runSamtoolsView(bamFilename, outFilename, bamOutput = True,
                    optionStr = "", verbose=False):
    """Runs samtools view, outputting to specified filename."""

    bamFlag = ""
    if bamOutput:
        bamFlag = "-b "
    cmdStr = "samtools view " + bamFlag + bamFilename + " " + optionStr + \
        " > " + outFilename
    runCommand(cmdStr, verbose)


def runSamtoolsIndex(currChromBamFilename, verbose = False):
    cmdStr = "samtools index " + currChromBamFilename
    runCommand(cmdStr, verbose)


def runIGVbatch(batchFilename, verbose = True):
    cmdStr = "java -Xmx20000m -jar /home/tom/Work/externalProgs/IGV-2.3.11/igv.jar -b " + batchFilename
    #cmdStr = "java -Xmx20000m -jar /mnt/oldHome/tom/Work/externalProgs/IGV-2.3.34/igv.jar -b " + batchFilename
    #cmdStr = "igv -b " + batchFilename
    runCommand(cmdStr, verbose)


def runBedtoolsIntersect(features1, features2, outputFilename, verbose = False):
    cmdStr = "bedtools intersect -a " + features1 + " -b " + features2 + " > " + outputFilename
    runCommand(cmdStr, verbose)


def runBigWigAverageOverBed(bedFilename, bigWigFilename, outputFilename, outputFilename2, verbose = False):
    cmdStr = "bigWigAverageOverBed -bedOut=" + outputFilename + " " + bigWigFilename + " " + bedFilename + " " + outputFilename2
    runCommand(cmdStr, verbose)


def runBigWigToWig(bigWigFilename, wigFilename, verbose=True):
    cmdStr = "bigWigToWig " + bigWigFilename + " " + wigFilename
    runCommand(cmdStr, verbose)


def runCommand(cmdStr, verbose, dieOnError = True):

    if verbose:
        print >> sys.stderr, "Running command:\n", cmdStr

    cmdResult = commands.getstatusoutput(cmdStr)

    if dieOnError and cmdResult[0] != 0:
        raise Exception("Command failed:\n" + cmdResult[1])

    if verbose:
        print >> sys.stderr, "Command result:", cmdResult


def run_MACS(chipBam, controlBam, outDir, name = "NoName", macsOptionsStr = "", verbose=False):
    """Runs my modified version of MACS."""

    startingDir = os.getcwd()

    # Figure out the absolute paths of the specified chip and control bam
    # filenames if they are not already absolute paths:
    if chipBam[0] == "/":
        absChipBamPath = chipBam
    else:
        absChipBamPath = startingDir + "/" + chipBam

    # If control is None, then use MACS2, with a hard-coded q-value threshold and no control
    # file specified (it's a quick hack):
    if controlBam == None:
        os.chdir(outDir)

        # Run the MACS command:
        cmdStr = "macs2 callpeak -t " + absChipBamPath + " --gsize=hs --format=BAM --name " + name + " --keep-dup=1 -q 0.01" + macsOptionsStr
        runCommand(cmdStr, verbose)

        # Change back to the former directory:
        os.chdir(startingDir)
        return

    if controlBam[0] == "/":
        absControlBamPath = controlBam
    else:
        absControlBamPath = startingDir + "/" + controlBam

    # Change to the specified directory:
    os.chdir(outDir)

    # Run the MACS command:
    cmdStr = "/mnt/oldHome/tom/Work/tomprogs/macs14.py -t " + absChipBamPath + " -c " + absControlBamPath + " --gsize=hs --format=BAM --mfold=2,6 --name " + name + " --keep-dup=1 " + macsOptionsStr
    runCommand(cmdStr, verbose)

    # Change back to the former directory:
    os.chdir(startingDir)
    

def deleteFiles(filenameList, verbose=False):
    """Deletes the specified list of files."""

    filenameListString = reduce(lambda x, y: x + " " + y, filenameList)
    cmdStr = "rm " + filenameListString
    runCommand(cmdStr, verbose)


def run_c_curve(bedFilename, outFilename, optionsStr="", verbose=False):
    """Runs the preseq c_curve program to estimate library complexity on the
    specified sorted input bed file."""

    cmdStr = "c_curve -o " + outFilename + " " + optionsStr + " " + bedFilename
    runCommand(cmdStr, verbose)


def run_lc_extrap(bedFilename, outFilename, optionsStr="", verbose=False):
    """Runs the preseq lc_extrap program to estimate library complexity on the
    specified sorted input bed file."""

    cmdStr = "lc_extrap -o " + outFilename + " " + optionsStr + " " + bedFilename
    runCommand(cmdStr, verbose)
