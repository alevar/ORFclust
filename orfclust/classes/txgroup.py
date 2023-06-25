from intervaltree import Interval, IntervalTree
from typing import Iterator

from classes.treader import TReader
from classes.transcript import Transcript, Object
from utils.common import *

class TXGroup:
    def __init__(self):
        self.objects = list()
        self.tid_map = dict() # transcript_id to position in objects list
    def clear(self):
        self.objects = list()
        self.tid_map = dict()
        
    def add_object(self,obj) -> int:
        idx = None
        if obj.get_type() in [Types.Transcript, Types.Exon, Types.CDS]:
            idx = self.tid_map.setdefault(obj.get_tid,len(self.objects))
            if len(self.objects) == idx:
                self.objects.append(Transcript(obj))
            self.objects[idx].add(obj)
        elif obj.get_type() in [Types.Bundle, Types.Gene]:
            for tx in obj.transcript_it():
                idx = self.tid_map.setdefault(tx.get_tid,len(self.objects))
                if len(self.objects) == idx:
                    self.objects.append(Transcript(tx))
                self.objects[idx].add(tx)
        else:
            raise Exception("Wrong object type provided to the add_object methods of TXGroup")
        
        return idx
    
    def __getitem__(self,idx):
        return self.objects[idx]
        
    def object_it(self):
        for obj in self.objects:
            yield obj

    def is_empty(self) -> bool:
        return len(self.objects) == 0

    def sort(self) -> None:
        # TODO: do we need other sorting methods?
        cmp = lambda obj: obj.get_cds()
        self._sort(cmp)
    
    def _sort(self,cmp) -> None:
        self.objects.sort(key=cmp)


    def to_gtf(self):
        res = ""
        for obj in self.objects:
            res+=obj.to_gtf()
        return res
    
    def to_gff(self):
        res = ""
        for obj in self.objects:
            res+=obj.to_gtf()
        return res
    
    __repr__ = to_gtf

    __str__ = __repr__
    

# largest specialization of the TXGroup which can yield all other types
class Transcriptome (TXGroup):
    def __init__(self):
        super().__init__()

    def build_from_file(self,fname: str) -> None:
        # get object (can e anything)
        # add to the transcriptome
        # transcriptome looks up by transcript id
        # what if object is a gene or a collection of objects?
        # neex to brake down into individual components with transcript IDs and add them individually
        treader = TReader(fname)
        for obj in treader.next_obj():
            self.add_object(obj)
        
    def coordinate_sort(self):
        self.objects.sort(key = lambda obj: obj.get_exons())
    def gid_sort(self):
        self.objects.sort(key = lambda obj: obj.get_gid())
            
    def transcript_it(self):
        for obj in self.object_it(): # TODO: types are enumeration shared between all classes
            if obj.type==Types.Transcript:
                yield obj
                
    def gene_it(self):
        gene = Gene()
        for obj in self.object_it():
            if gene.is_empty() or gene.get_gid() == obj.get_gid():
                gene.add_object(obj)
            else:
                yield gene
                gene = Gene()
                gene.add_object(obj)
    
        if not gene.is_empty():
            yield gene
            
    def bundle_it(self):
        bundle = Bundle()
        for obj in self.object_it():
            if bundle.is_empty() or bundle.overlaps(obj):
                bundle.add(obj)
            else:
                yield bundle
                bundle = Bundle()
                bundle.add(obj)
    
        if not bundle.is_empty():
            yield bundle

# specialization of a Object which guarantees that all objects overlap transitevily. Inherits also from the common object traits such as start,end,overlaps, etc
class Bundle (TXGroup,Object):
    def __init__(self):
        super().__init__(self)

# specialization of a Object which guarantees that all objects have the same gene ID. Inherits also from the common object traits such as start,end,overlaps, etc
class Gene (TXGroup,Object):
    def __init__(self):
        TXGroup.__init__(self)
        Object.__init__(self)
        
        self.intervals = IntervalTree()  # union of all exons in the locus (minus the introns)

    def add_object(self,obj: Object) -> bool:
        idx = TXGroup.add_object(self,obj)
        
        assert idx is not None,"wrong index in add_object of Gene"
        inserted_obj = self.objects[idx]

        self.intervals.update(inserted_obj.get_exons())
        self.intervals.merge_overlaps()


    # adds a transcript to the current locus
    def add_tx(self, tx: Transcript) -> None:
        assert self.seqid is None or self.seqid == tx.get_seqid(), "mismatching seqids: "+tx.get_tid()
        assert self.strand is None or self.strand == tx.get_strand(), "mismatching strands: "+tx.get_tid()

        self.seqid = tx.get_seqid()
        self.strand = tx.get_strand()

        self.intervals.update(tx.exons)
        self.intervals.merge_overlaps()

        self.start = min(self.start, tx.start)
        self.end = max(self.end, tx.end)

        self.txs.append(tx)

        return
    
    def get_start(self) -> int:
        return self.start
    def get_end(self) -> int:
        return self.end
    def get_gid(self) -> str:
        return self.objects[0].get_gid()