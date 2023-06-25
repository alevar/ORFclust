# TODO:
# general implementation of an independent annotation unit
# right now only works for transcript
# in the future it should work for constructing any type: transcript, gene, mRNA, etc - whatever an annotation may contain

# perhaps we can create a factory Object (like gffObj in gclib) which can create transcripts and other objects
# all objects will then be required to have some common methods that can be used by other operations (such as sorting). Methods would be: gene_id, start, end, etc

from intervaltree import Interval, IntervalTree
from utils.common import *


# every GTF object has some shared properties
# every object has coordinate range seqid, strand, start and end.
# every object can not span multiple chromosomes and seqids
# every object has source (otherwise set as default)
# every object can (but does not have to) have attributes
class Object:
    def __init__(self):
        self.start = float("inf") # inclusive
        self.end = 0 # non-inclusive

        self.source = "ORFclust"
        self.obj_type = None
        
        self.seqid = None
        self.strand = None

        self.attrs = dict()

    def clear(self):
        self.seqid = None
        self.strand = None
        self.source = "ORFclust"
        self.obj_type = None
        self.start = float("inf")
        self.end = 0
        self.attrs.clear()

    def is_empty(self) -> bool:
        return self.start == float("inf")

    def overlaps(self, obj) -> bool:
        return (self.get_start()<=obj.get_end() and self.get_end()>=obj.get_start())
    
    def set_seqid(self, seqid: str) -> None:
        self.seqid = seqid
    def set_strand(self, strand: str) -> None:
        self.strand = strand
    def set_source(self, source: str) -> None:
        self.source = source
    def set_type(self, obj_type: str) -> None:
        self.obj_type = obj_type
    def set_start(self, start: int):
        self.start = start
    def set_end(self, end: int):
        self.end = end
    def set_attributes(self, attrs: dict) -> None:
        self.attrs = attrs
    def add_attribute(self, k: str, v: str, replace: bool=False, append: bool=False) -> None: # replace - if key exists - replace it with new. append - if key exists - add new values to the list of values. If both enabled - will replace
        if append and not replace:
            if k in self.attrs:
                old_val = self.attrs[k]
                if not old_val == v:
                    self.attrs[k] = old_val+", "+v # TODO: should we replace with a list instead?
        else:
            self.attrs.setdefault(k,v)
            if replace:
                self.attrs[k] = v
    
    def get_source(self) -> str:
        return self.source
    def get_type(self) -> str:
        return self.obj_type
    def get_seqid(self) -> str:
        return self.seqid
    def get_strand(self) -> str:
        return self.strand
    def get_start(self) -> int:
        return self.start
    def get_end(self) -> int:
        return self.end
    def get_attr(self, attr: str) -> str:
        return self.attrs.get(attr ,None)
    def get_attributes(self) -> dict:
        return self.attrs
    
    # return none if any errors enountered
    # if store_all_attributes is True, then exon and CDS attributes will also be stored
    def add_line(self,gtf_line:str,store_all_attributes=False) -> bool:
        lcs = gtf_line.strip().split("\t")
        if not len(lcs) == 9:
            return None

        self.attrs = extract_attributes(lcs[8])

        lstart = int(lcs[3])
        lend = int(lcs[4])+1 # +1 here because intervaltree inclusive of the lower limit, but non-inclusive of the upper limit
        
        self.start = min(self.start,lstart)
        self.end = max(self.end,lend)
        self.seqid = lcs[0]
        self.strand = lcs[6]

        if lcs[2] == "transcript":
            self.obj_type = Types.Transcript
        elif lcs[2] == "exon":
            self.obj_type = Types.Exon
        elif lcs[2] == "CDS":
            self.obj_type = Types.CDS
        else:
            raise Exception("Wrong object type passed to Object.add_line()")
        
        return True
    
    def copy(self):
        obj = Object()
        obj.set_seqid(self.seqid)
        obj.set_source(self.source)
        obj.set_strand(self.strand)
        obj.set_attributes(self.attrs)
        return obj

class Transcript (Object):
    def __init__(self):
        super().__init__()

        self.obj_type = Types.Transcript
        self.exons = IntervalTree()
        self.cds = IntervalTree()

    def __init__(self, obj: Object):
        super().__init__()

        self.obj_type = Types.Transcript
        self.seqid = obj.seqid
        self.strand = obj.strand
        self.source = obj.source

        self.exons = IntervalTree()
        self.cds = IntervalTree()

        if obj.get_type() == Types.Transcript:
            self.from_transcript(obj)
        elif obj.get_type() == Types.Exon or obj.get_type() == Types.CDS:
            self.from_exon(obj)
        else:
            raise Exception("Wrong object type passed to Transcript constructor")

    def from_transcript(self, obj: Object) -> None:
        assert obj.get_type() == Types.Transcript,"wrong object type passed to Transcript constructor"
        
        self.attrs = rename_attributes(obj.attrs,{"ID":"transcript_id","Parent":"gene_id"})

        assert "transcript_id" in self.attrs,"transcript_id not found in attributes"
        assert "gene_id" in self.attrs,"gene_id not found in attributes"
        self.tid = self.attrs.get("transcript_id",None)
        self.gid = self.attrs.get("gene_id",None)

        self.start = obj.start
        self.end = obj.end

    def from_exon(self, obj: Object) -> None:
        assert obj.get_type() in [Types.Exon, Types.CDS],"wrong object type passed to Transcript constructor"
        
        self.attrs = rename_attributes(obj.attrs,{"Parent":"transcript_id"})

        assert "transcript_id" in self.attrs,"transcript_id not found in attributes"
        self.tid = self.attrs.get("transcript_id",None)
        if "gene_id" in self.attrs:
            self.gid = self.attrs.get("gene_id",None)

        self.start = obj.start
        self.end = obj.end

    def add(self,obj: Object) -> bool: # returns True if sucessfully added
        if self.strand != obj.get_strand() or self.seqid != obj.get_seqid():
            return False
        if self.tid != obj.get_attr("transcript_id"):
            return False
        
        if obj.get_type() == Types.Transcript:
            self.start = min(self.start,obj.get_start())
            self.end = max(self.end,obj.get_end())
            for k,v in obj.get_attributes().items():
                self.add_attribute(k,v)
        elif obj.get_type() == Types.Exon:
            self.start = min(self.start,obj.get_start())
            self.end = max(self.end,obj.get_end())
            self.exons.addi(obj.get_start(),obj.get_end(),obj)
        elif obj.get_type() == Types.CDS:
            self.cds.addi(obj.get_start(),obj.get_end(),obj)

        return True

    # makes sure that the transcript is complete
    # consider that only a transcript line has been added - this function will create required exons and intervals
    # if exon boundaries and transcript boundaries are out of sync (transcript line describes longer interval than exon chain) - they will be adjusted accordingly
    def finalize(self, obj: Object, extend: bool=False) -> None: # if extend is enabled - will extend 3' and 5' exons to match the start and end of transcrip, otherwise the transcript start end will be set to match existing start/end exon coordinates instead
        assert self.tid == self.attrs["transcript_id"],"transcript_id missing - can not assign"
        if len(self.exons)>0:
            exon_start = sorted(self.exons)[0][0]
            exon_end = sorted(self.exons)[-1][1]
            if extend:
                self.start = min(self.start,exon_start)
                self.end = max(self.end,exon_end)
                self.exons.smallest = self.start # TODO: implement
                self.exons.largest = self.end # TODO: implement
            else:
                self.start = exon_start
                self.end = exon_end
        else:
            obj = Object()
            obj.set_seqid(self.seqid)
            obj.set_strand(self.strand)
            obj.set_start(self.start)
            obj.set_end(self.end)
            obj.set_type("exon")
            obj.set_attributes({"transcript_id":self.tid})
            exon = Exon(obj)
            self.exons.addi(self.start,self.end,exon)

        # make sure the CDS (if added) fits within exons
        if len(self.cds)>0:
            assert len(self.exons[self.cds.begin()])>0 and len(self.exons[self.cds.end()])>0,"invalid CDS in detected in transcript when finalizing: "+self.tid

    def clear(self):
        super().clear()

        self.exons = IntervalTree()
        self.cds = IntervalTree()

        self.tid = None

    def set_exons(self ,exons: list[tuple[int,int]]) -> None:
        self.exons = IntervalTree.from_tuples(exons)
    def set_cds(self ,cds: list[tuple[int,int]]) -> None:
        self.cds = IntervalTree.from_tuples(cds)
    def set_tid(self ,tid:str) -> None:
        self.tid = tid
    def set_gid(self,gid:str) -> None:
        self.gid = gid

    def nume(self):
        return len(self.exons)
    def numc(self):
        return len(self.cds)
    def has_cds(self):
        return len(self.cds) > 0
    # yields all introns of the transcript
    # intervals are also inclusive of the lowe bound but non-inclusive of the upper bound
    # [first position of the intron, last position of the intron)
    def get_introns(self):
        if len(self.exons) > 1:
            prev_exon = None
            for e in sorted(self.exons):
                if prev_exon is not None:
                    yield Interval(prev_exon[1], e[0])
                prev_exon = e
    def get_exons(self):
        return sorted(self.exons)
    def get_cds(self):
        return sorted(self.cds)
    def get_tid(self):
        return self.tid
    def get_gid(self):
        return self.gid
    def get_cstart(self):
        return sorted(self.cds)[0][0]
    def get_cend(self):
        return sorted(self.cds)[-1][1]

    def to_gtf(self):
        res =  self.seqid+"\t"+\
               self.source+"\t"+\
               Types.type2str(self.obj_type) +"\t"+\
               str(self.start)+"\t"+\
               str(self.end-1)+"\t"+\
               "."+"\t"+\
               self.strand+"\t"+\
               "."+"\t"+\
               to_attribute_string(self.attrs,False,"transcript")+"\n"
        
        for e in sorted(self.exons):
            res += e.to_gtf() + "\n"
            
        for c in sorted(self.cds):
            res += c.to_gtf() + "\n"

        return res.rstrip("\n")

    def to_gff(self):
        res =  self.seqid+"\t"+\
               self.source+"\t"+\
               Types.type2str(self.obj_type) +"\t"+\
               str(self.start)+"\t"+\
               str(self.end-1)+"\t"+\
               "."+"\t"+\
               self.strand+"\t"+\
               "."+"\t"+\
               to_attribute_string(self.attrs,True,"transcript")+"\n"
        
        for e in sorted(self.exons):
            res += e.to_gff() + "\n"
            
        for c in sorted(self.cds):
            res += c.to_gff() + "\n"

        return res.rstrip("\n")

    __repr__ = to_gtf

    __str__ = __repr__

    def copy(self):
        tx = Transcript()
        tx.set_tid(self.tid)
        tx.set_gid(self.gid)
        tx.set_exons(self.exons)
        tx.set_cds(self.cds)
        return tx
    
class Exon(Object):
    def __init__(self):
        super().__init__()
        self.tid = None
        self.gid = None

        self.obj_type = Types.Exon

    def __init__(self,obj: Object):
        super().__init__()
        self.tid = None
        self.gid = None

        self.obj_type = Types.Exon
        self.seqid = obj.seqid
        self.strand = obj.strand
        self.source = obj.source
        
        if obj.get_type() == Types.Exon:
            self.from_exon(obj)
        else:
            raise Exception("wrong object type passed to Exon constructor")

    def from_exon(self, obj):
        assert obj.get_type() == Types.Exon or obj.get_type() ==  Types.CDS,"wrong object type passed to Transcript constructor"
        
        self.attrs = rename_attributes(obj.attrs,{"Parent":"transcript_id"})

        assert "transcript_id" in self.attrs,"transcript_id not found in attributes"
        self.tid = self.attrs.get("transcript_id",None)
        if "gene_id" in self.attrs:
            self.gid = self.attrs.get("gene_id",None)

        self.start = obj.start
        self.end = obj.end

    def set_tid(self ,tid:str) -> None:
        self.tid = tid
    def set_gid(self,gid:str) -> None:
        self.gid = gid

    def get_tid(self):
        return self.tid
    def get_gid(self):
        return self.gid
    def get_cstart(self):
        return sorted(self.cds)[0][0]
    def get_cend(self):
        return sorted(self.cds)[-1][1]

    def to_gtf(self):
        res = self.seqid+"\t"+\
                self.source+"\t"+\
                Types.type2str(self.obj_type) +"\t"+\
                str(self.start)+"\t"+\
                str(self.end-1)+"\t"+\
                "."+"\t"+\
                self.strand+"\t"+\
                "."+"\t"+\
                "transcript_id \""+self.tid+"\";"
        if self.gid is not None:
            res += " gene_id \""+self.gid+"\";"
        
        res += to_attribute_string(self.attrs,False,"exon")
        return res
    
    def to_gff(self):
        res = self.seqid+"\t"+\
                self.source+"\t"+\
                Types.type2str(self.obj_type) +"\t"+\
                str(self.start)+"\t"+\
                str(self.end-1)+"\t"+\
                "."+"\t"+\
                self.strand+"\t"+\
                "."+"\t"+\
                "Parent="+self.tid+";"
        if self.gid is not None:
            res += " gene_id="+self.gid+";"
        
        res += to_attribute_string(self.attrs,True,"exon")
        return res
    
    __repr__ = to_gtf

    __str__ = __repr__

class CDS(Exon):
    def __init__(self):
        Object.__init__(self)
        self.tid = None
        self.gid = None

        self.obj_type = Types.CDS

    def __init__(self,obj: Object):
        Object.__init__(self)
        self.tid = None
        self.gid = None

        self.obj_type = Types.CDS
        self.seqid = obj.seqid
        self.strand = obj.strand
        self.source = obj.source
        
        if obj.get_type() == Types.CDS:
            self.from_exon(obj)
        else:
            raise Exception("wrong object type passed to CDS constructor")

class GTFObjectFactory:
    @staticmethod
    def create(line):
        obj = Object()
        obj.add_line(line)
        if obj.get_type() == Types.Transcript:
            return Transcript(obj)
        elif obj.get_type() == Types.Exon:
            return Exon(obj)
        elif obj.get_type() == Types.CDS:
            return CDS(obj)
        else:
            return obj