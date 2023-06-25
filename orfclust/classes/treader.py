from typing import Iterator

from classes.transcript import GTFObjectFactory, Object
from utils.common import *

class TReader:
    def __init__(self,fname:str=None):
        self.fname = ""
        self.fp = None
        self.gff = True # set to gtf if the file format is detected to be GTF, otherwise is set ot False for GFF
        self.comments = [] # contains a list of all comment lines encountered in the gtf file (including leading #)
        
        if not fname is None:
            self.set_file(fname)

    def __del__(self):
        if not self.fp is None:
            self.fp.close()

    def set_file(self,fname:str):
        self.fname = fname
        if self.fp is not None:
            self.fp.close()

        self.gff = self.is_gtf()
        if self.gff is None:
            raise Exception("File is not a valid GTF or GFF file: "+fname)
        self.fp = open(self.fname,"r")

    # function opens the file, peaks inside, and returns True if the file is in GTF format, False if it is in GFF format and None if it is not a valid file
    # will always close the file and re-open it
    def is_gtf(self) -> bool:
        gff = None
        if not self.fp is None:
            self.fp.close()

        self.fp = open(self.fname,"r")

        for line in self.fp:
            if line.startswith("#"):
                continue
            lcs = line.strip().split("\t")
            if len(lcs) > 9:
                gff = None
                break
            
            if lcs[2] not in ["transcript","exon","CDS"]:
                gff = None
                break

            if lcs[2] == "transcript":
                if lcs[8].startswith("ID="):
                    gff = True
                    break
                elif lcs[8].startswith("transcript_id"):
                    gff = False
                    break
                else:
                    gff = None
                    break
            else:
                continue
        
        self.fp.close()
        self.fp = open(self.fname,"r")
        return gff

    def next_obj(self) -> Iterator[Object]:
        for line in self.fp:
            if line.startswith('#'):
                self.comments.append(line)
                continue

            obj = GTFObjectFactory.create(line)
            if not obj is None:
                yield obj

# needs a next_group method
# the next_group method should yield groups of transcripts based on some metric. The metric should be specified by a comparator function passed to the method
# for example, we can have a comparator function which returns true if two transcripts have the same ORF and returns false otherwise


# transcriptome should be responsible for reading the data
# perhaps rename into treader instead
# the goal is to load transcripts into loci (whether by


# since we are working with loci instead of a full transcriptome, perhaps isntead of the "yield or next_group" method, we can simply have a method that returns all groups as a collection
# we can then choose how to iterate or what to do with those groups separately


# then the class structure will be
# TReader
# Locus
# Transcript