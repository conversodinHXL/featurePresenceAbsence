#!/usr/bin/env python

import os
import sys
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIXML
import threading
import argparse
import subprocess
import time

programUsage = """
NAME
	preAb.py - Presence/absence of genes with BLAST.
SYNOPSIS
	preAb.py -i genes*.fasta -d db*.fasta -o output -c query_cov -p query_identity
	-e E-value
DESCRIPTION
	Runs BLAST against the query genes against the given database sequences to identify
	whether it's present/absent.

	By default the format for the input sequences is FASTA.

	#Simplest way to run it is to provide the input sequences in fasta format (default).
	pneumoSerotyper.py -i *.fasta -o output
OPTIONS
	-h
	Help. Print this help information and exit.

	-i filenames
	Specify the input fasta sequence files for the program. You can use wildcard (*)
	to speficy multiple files.

	-o filename
	Specify the prefix for the output files. By default the output files will have
	a prefix "Serotypes.Summary" if not specified by the user.

	-c value
	Specify the minimum query coverage (%) to be used.

	-p value
	Specify the minimum identity (%) for all matched regions
AUTHOR
	Chrispin Chaguza, Chrispin.Chaguza@liverpool.ac.uk. June 2015

DEPENDENCIES
	biopython,NCBI BLAST
"""


def readUserArguments(UserArgs):
    Options = argparse.ArgumentParser(UserArgs[0],
                    description="preAb.py - Presence/absence of genes with BLAST",
                    prefix_chars="-",
                    add_help=False,
                    epilog="Chrispin Chaguza (Chrispin.Chaguza@liv.ac.uk)")

    Options.add_argument("-i",
                    action="store",
                    nargs="*",
                    required=False,
                    metavar="Input",
                    dest="Input",
                    help="Input genomes (fasta format)")

    Options.add_argument("-d",
                    action="store",
                    nargs="*",
                    required=False,
                    metavar="Dbase",
                    dest="Dbase",
                    help="Input genomes (fasta format)")


    Options.add_argument("-o",
                     action="store",
                     nargs=1,
                     required=False,
                     metavar="Output",
                     dest="Output",
                     help="Output prefix (default=Output)",
                     default="Output")

    Options.add_argument("-c",
                     action="store",
                     nargs=1,
                     required=False,
                     metavar="Sequence_Coverage",
                     dest="Sequence_Coverage",
                     help="Percent query coverage (default=90)",
                     default="90")

    Options.add_argument("-p",
                     action="store",
                     nargs=1,
                     required=False,
                     metavar="Sequence_Identity",
                     dest="Sequence_Identity",
                     help="Percent Identity for BLAST HSPs (default=90)",
                     default="90")

    Options.add_argument("-r",
                     action="store_true",
                     dest="Keep",
                     help="Remove BLAST output files (keep by default)")

    Options.add_argument("-h",
                     action="store_true",
                     dest="Help",
                     help="Show detailed help")

    Options = Options.parse_args()

    return Options


def checkUserArguments(UserOptions):
    Options = UserOptions
    OptionsVars = {}

    if Options.Help:
        sys.stdout.write(str(programUsage) + "\n")
        sys.exit()

    if Options.Input:
        OptionsVars["i"] = Options.Input[0:]
    else:
        showErrorMessage("input files (-i) are required")
        sys.exit()

    if Options.Dbase:
        OptionsVars["d"] = Options.Dbase[0:]
    else:
        showErrorMessage("input files (-d) are required")
        sys.exit()

    if Options.Output != "Output":
        OptionsVars["o"] = Options.Output[0:][0]
    else:
        OptionsVars["o"] = Options.Output[0:]

    if Options.Sequence_Coverage != "90":
        if ((int(Options.Sequence_Coverage[0:][0]) > 0) and (int(Options.Sequence_Coverage[0:][0]) < 100)):
            OptionsVars["c"] = Options.Sequence_Coverage[0:][0]
        else:
            showErrorMessage("Sequence coverage (-c) should be >0 and <=100 (default=90)")
    else:
        OptionsVars["c"] = Options.Sequence_Coverage[0:]

    if Options.Sequence_Identity != "90":
        if ((float(Options.Sequence_Identity[0:][0]) > 0) and (float(Options.Sequence_Identity[0:][0]) <= 100)):
            OptionsVars["p"] = Options.Sequence_Identity[0:][0]
        else:
            showErrorMessage("Sequence identity (-p) should be >0 and <=100 (default=90)")
    else:
        OptionsVars["p"] = Options.Sequence_Identity[0:]

    OptionsVars["r"] = Options.Keep

    return OptionsVars


def showErrorMessage(ErrorMessage):
    sys.stdout.write("\nError: " + str(ErrorMessage) + "\n")
    sys.stdout.write("\nUse -h option to see more detailed help\n")

    sys.exit()


def showProgramStatus(ItemList, ItemPos):
    NumElements = ItemList
    ProgressChars = "="

    HashChars = ProgressChars * int(round(float(ItemPos + 1) / float(NumElements) * 100))
    SpaceChars = " " * int(round(100 - len(HashChars)))
    PercentChars = float(ItemPos + 1) / NumElements * 100

    if (ItemPos + 1) >= ItemList or PercentChars >= 100 or \
            (int(float(ItemPos + 1) / float(NumElements) * 100) >= 100):
        sys.stdout.write("\r|{0}| {1:.2f}% {2}".format(ProgressChars * 100, 100, ""))
        sys.stdout.flush()
        sys.stdout.write("\r|{0}| {1:.2f}% ".format(ProgressChars * 100, 100))

    else:
        sys.stdout.write("\r|{0}| {1:.2f}%  {2}".format(HashChars + SpaceChars, PercentChars, ""))
        sys.stdout.flush()
        sys.stdout.write("\r|{0}| {1:.2f}%  {2}".format(HashChars + SpaceChars, PercentChars, ""))


def RunBLASTThread(blastCommand, blastOutput):
    try:
        FNULL = open(blastOutput, "wb")
        subprocess.check_call(blastCommand, stdout=FNULL, stderr=FNULL)
        FNULL.close()

    except (StandardError, KeyboardInterrupt, SystemExit), ErrorText:
        sys.stdout.write("\n" + ErrorText.message + "\n\n")
        sys.exit()


def main():
    Args = readUserArguments(sys.argv[:])
    inputOptions = checkUserArguments(Args)

    print "\n---------------------------------------------------------------------------------------------------------------"
    print "-   preAb v1.0.0                                                                                              -"
    print "---------------------------------------------------------------------------------------------------------------"
    print "-   Program for identifying presence/absence of genes sequences using BLAST                                   -"
    print "-   Institute of Infection and Global Health, University of Liverpool, UK                                     -"
    print "-   All rights reserved                                                                                       -"
    print "---------------------------------------------------------------------------------------------------------------"


    if not os.path.exists(inputOptions["o"] + ".Tmp.Files/"):
        os.makedirs(inputOptions["o"] + ".Tmp.Files/")
    else:
        pass

    print "\ncreating sequence database for the input fasta sequences"

    for geneSeqPos, geneSeq in enumerate(inputOptions["d"]):
        makeDBaseCommand = []
        makeDBaseCommand.append("makeblastdb")
        makeDBaseCommand.append("-in")
        makeDBaseCommand.append(geneSeq)
        makeDBaseCommand.append("-dbtype")
        makeDBaseCommand.append("nucl")

        FNULL = open(os.devnull, "wb")
        subprocess.call(makeDBaseCommand, stdout=FNULL, stderr=FNULL)
        FNULL.close()

        showProgramStatus(len(inputOptions["d"]), geneSeqPos)

    print "\n\nrunning BLAST on the query sequences against database"

    blastOutputFiles = []
    inputGeneNames = []
    inputDBNames = []

    currentThreads = threading.activeCount()

    for geneSeqPos, geneSeq in enumerate(inputOptions["d"]):
        inSeqFName = os.path.basename(geneSeq).split(".")[0]
        inputDBNames.append(inSeqFName)

        for isolateSeqPos, isolateSeq in enumerate(inputOptions["i"]):
            geneFName = os.path.basename(isolateSeq).split(".")[0]

            blastCommand = []

            blastCommand.append("blastn")
            blastCommand.append("-query")
            blastCommand.append(geneSeq)
            blastCommand.append("-db")
            blastCommand.append(isolateSeq)
            blastCommand.append("-outfmt")
            blastCommand.append("5")
            blastCommand.append("-out")
            blastCommand.append(inputOptions["o"] + ".Tmp.Files/" + inSeqFName + "." + geneFName + ".blast.xml")

            blastOutputFiles.append(inputOptions["o"] + ".Tmp.Files/" + inSeqFName + "." + geneFName + ".blast.xml")
            inputGeneNames.append(geneFName)

            while True:
                if threading.activeCount() < 15:
                    try:
                        blastThread = threading.Thread(name=inSeqFName + "." + geneFName,
                                                       target=RunBLASTThread, args=(blastCommand, os.devnull))
                        blastThread.setDaemon(True)
                        blastThread.start()

                    except (StandardError, KeyboardInterrupt, SystemExit):
                        print "Unknown error occurred OR the user killed the program"
                        sys.exit()

                    showProgramStatus(len(inputOptions["i"])*len(inputOptions["d"]), geneSeqPos*len(inputOptions["i"])+isolateSeqPos)

                    break

        showProgramStatus(len(inputOptions["i"])*len(inputOptions["d"]), (geneSeqPos+1)*len(inputOptions["i"]))

    while threading.activeCount() > currentThreads:
        time.sleep(2)

    print "\n\ncalculating length of each query sequence"

    blastOutputTXT = open(inputOptions["o"]+".Summary.txt","w")

    geneSeqLengths = {}

    for geneSequencePos,geneSequence in enumerate(inputOptions["i"]):
        geneSequenceRecord=SeqIO.read(open(geneSequence,"rU"),"fasta")
        geneSeqLengths[os.path.basename(geneSequence).split(".")[0]]=len(geneSequenceRecord.seq)

        showProgramStatus(len(inputOptions["i"]), geneSequencePos)

    blastOutputTXT.write("taxa")

    for querySeqPos, querySeq in enumerate(inputOptions["i"]):
        blastOutputTXT.write("\t"+os.path.basename(querySeq).split(".")[0])

    blastOutputTXT.write("\n")

    print "\n\ndetermining Cps sequence coverage and nucleotide identity"

    for subjectSeqPos, subjectSeq in enumerate(inputDBNames):

        blastOutputTXT.write(subjectSeq)
        sbjctSeqName=os.path.basename(subjectSeq).split(".")[0]

        for querySeqPos, querySeq in enumerate(inputOptions["i"]):

            querySeqName=os.path.basename(querySeq).split(".")[0]

            try:

                blastOutput=open(inputOptions["o"] + ".Tmp.Files/"+sbjctSeqName+"."+querySeqName+".blast.xml","rU")
                blastRecord=NCBIXML.read(blastOutput)

                tempCoverage = []
                tempIdentity = []
                tempIdentityCount = 0

                for blastRecordAlign in blastRecord.alignments:
                    for hsp in blastRecordAlign.hsps:
                        if (hsp.expect < 0.0001) and (len(hsp.query) > 400):
                            tempCoverage.append(len(hsp.query))

                            tempIdentity.append(hsp.match.count("|")/float(len(hsp.query)))
                            tempIdentityCount+=1

                querySeqCoverage = (sum(tempCoverage)/float(geneSeqLengths[querySeqName]))*100

                if tempIdentityCount == 0:
                    matchIdentity = 0
                else:
                    matchIdentity = (sum(tempIdentity)/float(tempIdentityCount))*100

                if querySeqCoverage >= 100:
                    querySeqCoverage = 100
                else:
                    querySeqCoverage=querySeqCoverage

                if (querySeqCoverage>=float(inputOptions["c"])) and (matchIdentity>=float(inputOptions["p"])):
                    blastOutputTXT.write("\t1")
                else:
                    blastOutputTXT.write("\t0")

            except (KeyboardInterrupt, SystemExit),err:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname, exc_tb.tb_lineno)
                print "Something wrong with BLAST files OR the user killed the program"
                sys.exit()

            showProgramStatus(len(inputDBNames)*len(inputOptions["i"]),subjectSeqPos*len(inputOptions["i"])+querySeqPos)

        showProgramStatus(len(inputDBNames)*len(inputOptions["i"]),(subjectSeqPos+1)*len(inputOptions["i"]))

        #sortedSerotypeCoverage = sorted(serotypeCoverage.iteritems(),key=operator.itemgetter(1),reverse=True)
        blastOutputTXT.write("\n")

    blastOutputTXT.close()

    print "\n--------------------------------------------------------------------------------------------------------------"
    print "-        Finished                                                                                            -"
    print "--------------------------------------------------------------------------------------------------------------\n"

if __name__ == "__main__":
    main()
