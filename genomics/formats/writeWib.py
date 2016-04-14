#!/usr/bin/env python
# DESCRIPTION: A command line interface for writing wib files (converted
# from wig).

import commands, math, os, pdb, random, shutil, string, struct, sys
import bed, wibIO
from optparse import OptionParser

def main():
    # Parse the command-line arguments...
    description = """usage: %prog [options] <wigFile> <outDir>\n

Writes binary wiggle (wib) data to the specified directory, <wibDir>.
Reads data from <wigFile> specifying the wiggle data for a single
chromosome only. Infers the name of the chromosome (and hence the name of the
output file) from the chromosome name in the input wiggle file.
"""
#    pdb.set_trace()

    parser = OptionParser(usage = description)
    parser.add_option("--lower", dest = "lower", default = "0",
                      help = "The lower value of the range of wig values.")
    parser.add_option("--upper", dest = "upper", default = "1",
                      help = "The upper value of the range of wig values.")
    parser.add_option("--null", dest = "null", default = "-1",
                      help = "The null value of the wig values.")
    parser.add_option("--debug", action="store_true", dest="debug",
                      help = "Debug the program using pdb.")
    (options, args) = parser.parse_args()

    if (options.debug):
        pdb.set_trace()

    if (len(args) != 2):
        print >> sys.stderr, "Invalid number of input args:", len(args)
        parser.print_help()
        sys.exit(1)

    # Parse the input parameters...

    # Make a new converter from the specified range values:
    lower = float(options.lower)
    upper = float(options.upper)
    null = float(options.null)
    conv = wibIO.wibConverter(lower, upper, null)

    wigFile = args[0]
    outDir = args[1]

    # Set up a writer to write the wib data:
    wibRdr = wibIO.wibWriter(outDir, open(wigFile), conv)

    # Write the binary wib file:
    wibRdr.writeWib()

    # Write out the range file specifying the lower, upper, and null values:
    rangeFile = open(outDir + "/range.txt", 'w')
    print >> rangeFile, lower
    print >> rangeFile, upper
    print >> rangeFile, null


if __name__ == "__main__":
    sys.exit(main())
