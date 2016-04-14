#!/usr/bin/env python
# DESCRIPTION: A command line interface for reading wib files (and converting
# to wig).

import commands, math, os, pdb, random, shutil, string, struct, sys
import bed, wibIO
from optparse import OptionParser

def main():
    # Parse the command-line arguments...
    description = """usage: %prog [options] <wibDir>\n

Reads binary wiggle (wib) data from the specified directory, <wibDir>.
Retrieves data for each of the bed regions specified on stdin. Prints the
retrieved wiggle data to stdout, with each region's data in the following
format:

# Comment line, stating chromosome location
[data lines, one line per base-pair]+
"""

    parser = OptionParser(usage = description)
    parser.add_option("--considerStrand", dest = "considerStrand",
                      default = False, action = "store_true",
                      help = "Reverse direction of wiggle data for bed " + \
                          "regions on minus strand.")
    parser.add_option("--i", dest = "inFile",
                      default = None,
                      help = "Take input bed from specified file instead " + \
                          "of stdin.")
    parser.add_option("--debug", action="store_true", dest="debug",
                      help = "Debug the program using pdb.")
    # parser.add_option("--debug", action="store_true", dest="debug",
    #                   help = "Debug the program using pdb.")
    (options, args) = parser.parse_args()

    if (options.debug):
        pdb.set_trace()

    # if (options.debug):
    #     debugger = pdb.Pdb(stdin=sys.stdin)
    #     debugger.use_rawinput = False
    #     print >> sys.stderr, debugger
    #     debugger.set_trace()

    if (len(args) != 1):
        print >> sys.stderr, "Invalid number of input args:", len(args)
        parser.print_help()
        sys.exit(1)

    # Read bed file input information from stdin:
    if (options.inFile == None):
        inBEDs = bed.BED_Track(sys.stdin)
    else:
        inBEDs = bed.BED_Track(open(options.inFile))

    # Set up a reader to read the wib data:
    wibRdr = wibIO.wibReader(args[0])

    # Read the requested regions:
    wigData = wibRdr.getWigs(inBEDs, flipOnStrand=options.considerStrand)

    # Print them to stdout:
    print >> sys.stderr, "Printing wig data to stdout..."
    for bedLine in inBEDs.getAllBEDs():
        print >> sys.stdout, "#", bedLine.to_string()
        if (wigData[bedLine] != None):
            for wigVal in wigData[bedLine]:
                print >> sys.stdout, wigVal
    print >> sys.stderr, "Finished printing."


if __name__ == "__main__":
    sys.exit(main())
