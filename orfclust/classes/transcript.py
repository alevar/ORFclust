from intervaltree import Interval, IntervalTree
from utils.common import *

class Object:
    def __init__(self):
        self.start = float("inf") # inclusive
        self.end = 0 # non-inclusive

        self.source = "ORFclust"
        self.obj_type = None
        
        self.seqid = None
        self.strand = None

        self.attrs = dict()
        self.expression = list()

    def clear(self):
        """
        Clear the attributes of the object.
        """
        self.seqid = None
        self.strand = None
        self.source = "ORFclust"
        self.obj_type = None
        self.start = float("inf")
        self.end = 0
        self.attrs.clear()
        self.expression = list()

    def is_empty(self) -> bool:
        """
        Check if the object is empty.

        Returns:
            bool: True if the object is empty, False otherwise.
        """
        return self.start == float("inf")

    def overlaps(self, obj) -> bool:
        """
        Check if the object overlaps with another object.

        Args:
            obj (Object): Another object to check for overlap.

        Returns:
            bool: True if the objects overlap, False otherwise.
        """
        return (self.get_start()<=obj.get_end() and self.get_end()>=obj.get_start())
    
    def set_seqid(self, seqid: str) -> None:
        """
        Set the sequence ID of the object.

        Args:
            seqid (str): The sequence ID.

        Returns:
            None
        """
        self.seqid = seqid

    def set_strand(self, strand: str) -> None:
        """
        Set the strand of the object.

        Args:
            strand (str): The strand.

        Returns:
            None
        """
        self.strand = strand

    def set_source(self, source: str) -> None:
        """
        Set the source of the object.

        Args:
            source (str): The source.

        Returns:
            None
        """
        self.source = source

    def set_type(self, obj_type: str) -> None:
        """
        Set the type of the object.

        Args:
            obj_type (str): The object type.

        Returns:
            None
        """
        self.obj_type = obj_type

    def set_start(self, start: int):
        """
        Set the start position of the object.

        Args:
            start (int): The start position.

        Returns:
            None
        """
        self.start = start

    def set_end(self, end: int):
        """
        Set the end position of the object.

        Args:
            end (int): The end position.

        Returns:
            None
        """
        self.end = end

    def set_attributes(self, attrs: dict) -> None:
        """
        Set the attributes of the object.

        Args:
            attrs (dict): A dictionary of attributes.

        Returns:
            None
        """
        self.attrs = attrs

    def add_attribute(self, k: str, v: str, replace: bool=False, append: bool=False) -> None: # replace - if key exists - replace it with new. append - if key exists - add new values to the list of values. If both enabled - will replace
        """
        Add or update an attribute to the object.

        Args:
            k (str): The attribute key.
            v (str): The attribute value.
            replace (bool, optional): If True and the key exists, replace the value. Defaults to False.
            append (bool, optional): If True and the key exists, append the value to the list. Defaults to False.

        Returns:
            None
        """
        if append and not replace:
            if k in self.attrs:
                old_val = self.attrs[k]
                if not old_val == v:
                    self.attrs[k] = old_val+", "+v # TODO: should we replace with a list instead?
        else:
            self.attrs.setdefault(k,v)
            if replace:
                self.attrs[k] = v

    def add_expression(self, exp: float) -> None:
        """
        Add an expression value to the object.

        Args:
            exp (float): The expression value.

        Returns:
            None
        """
        self.expression.append(exp)

    def get_source(self) -> str:
        """
        Get the source of the object.

        Returns:
            str: The source.
        """
        return self.source
    
    def get_type(self) -> str:
        """
        Get the type of the object.

        Returns:
            str: The object type.
        """
        return self.obj_type

    def get_seqid(self) -> str:
        """
        Get the sequence ID of the object.

        Returns:
            str: The sequence ID.
        """
        return self.seqid

    def get_strand(self) -> str:
        """
        Get the strand of the object.

        Returns:
            str: The strand.
        """
        return self.strand

    def get_start(self) -> int:
        """
        Get the start position of the object.

        Returns:
            int: The start position.
        """
        return self.start

    def get_end(self) -> int:
        """
        Get the end position of the object.

        Returns:
            int: The end position.
        """
        return self.end

    def get_attr(self, attr: str) -> str:
        """
        Get the value of a specific attribute.

        Args:
            attr (str): The attribute name.

        Returns:
            str: The attribute value, or None if the attribute does not exist.
        """
        return self.attrs.get(attr, None)

    def get_attributes(self) -> dict:
        """
        Get all attributes of the object.

        Returns:
            dict: A dictionary of attributes.
        """
        return self.attrs

    def get_expression(self) -> list:
        """
        Get the expression values of the object.

        Returns:
            list: A list of expression values.
        """
        return self.expression

    def add_line(self,gtf_line:str,store_all_attributes=False) -> bool:
        """
        Add information from a GTF line to the object.

        Args:
            gtf_line (str): The GTF line.
            store_all_attributes (bool, optional): If True, store all attributes including exon and CDS attributes.
                                                  Defaults to False.

        Returns:
            bool: True if the line was successfully added, False otherwise. Returns None in case of errors
        """
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
        """
        Create a copy of the object.

        Returns:
            Object: A new instance of the Object class with the same attribute values.
        """
        obj = Object()
        obj.set_seqid(self.seqid)
        obj.set_source(self.source)
        obj.set_strand(self.strand)
        obj.set_attributes(self.attrs)
        obj.set_type(self.obj_type)
        obj.set_attributes(self.attrs)
        obj.set_start(self.start)
        obj.set_end(self.end)
        obj.expression = [x for x in self.expression]
        return obj

class Transcript (Object):
    """
    A class representing a transcript.

    Attributes:
        obj_type (Types): The type of the object (Transcript).
        exons (IntervalTree): An interval tree containing exons.
        cds (IntervalTree): An interval tree containing CDS regions.

    Inherits from:
        Object

    """
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
        """
        Initialize the Transcript object from a general object with Transcript object.

        Args:
            obj (Object): Object with type Transcript.

        Returns:
            None

        """
        assert obj.get_type() == Types.Transcript,"wrong object type passed to Transcript constructor"
        
        self.attrs = rename_attributes(obj.attrs,{"ID":"transcript_id","Parent":"gene_id"})

        assert "transcript_id" in self.attrs,"transcript_id not found in attributes"
        assert "gene_id" in self.attrs,"gene_id not found in attributes"
        self.tid = self.attrs.get("transcript_id",None)
        self.gid = self.attrs.get("gene_id",None)

        self.start = obj.start
        self.end = obj.end

    def from_exon(self, obj: Object) -> None:
        """
        Initialize the Transcript object from an Exon or CDS object.

        Args:
            obj (Object): The Exon or CDS object.

        Returns:
            None

        """
        assert obj.get_type() in [Types.Exon, Types.CDS],"wrong object type passed to Transcript constructor"
        
        self.attrs = rename_attributes(obj.attrs,{"Parent":"transcript_id"})

        assert "transcript_id" in self.attrs,"transcript_id not found in attributes"
        self.tid = self.attrs.get("transcript_id",None)
        if "gene_id" in self.attrs:
            self.gid = self.attrs.get("gene_id",None)

        self.start = obj.start
        self.end = obj.end

    def add(self,obj: Object) -> bool: # returns True if sucessfully added
        """
        Add an Object to the Transcript. Objects are sorted accordingly to update coordinates, boundaries, exon and cds chains

        Args:
            obj (Object): The Object to add.

        Returns:
            bool: True if the Object was successfully added, False otherwise.

        """
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

    def finalize(self, extend: bool=False) -> None:
        """
        Make sure that the transcript is complete.

        If extend is enabled, the 3' and 5' exons will be extended to match the start and end of the transcript.
        # consider that only a transcript line has been added - this function will create required exons and intervals
        # if exon boundaries and transcript boundaries are out of sync (transcript line describes longer interval than exon chain) - they will be adjusted accordingly

        Args:
            extend (bool, optional): Whether to extend the exons. Defaults to False.

        Returns:
            None

        """
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
            obj.set_type(Types.Exon)
            obj.set_attributes({"transcript_id":self.tid})
            exon = Exon(obj)
            self.exons.addi(self.start,self.end,exon)

        # make sure the CDS (if added) fits within exons
        if len(self.cds)>0:
            assert len(self.exons[self.cds.begin()])>0 and len(self.exons[self.cds.end()-1])>0,"invalid CDS detected in transcript when finalizing: "+self.tid

    def clear(self):
        """
        Clear the Transcript object.

        Returns:
            None

        """
        super().clear()

        self.exons = IntervalTree()
        self.cds = IntervalTree()
        self.tid = None
        self.gid = None
        self.obj_type = Types.Transcript

    def set_exons(self, exons: list[tuple[int, int]]) -> None:
        """
        Set the exons of the Transcript.

        Args:
            exons (list[tuple[int, int]]): A list of exon intervals represented as tuples.

        Returns:
            None

        """
        self.exons = IntervalTree.from_tuples(exons)

    def set_cds(self, cds: list[tuple[int, int]]) -> None:
        """
        Set the CDS regions of the Transcript.

        Args:
            cds (list[tuple[int, int]]): A list of CDS intervals represented as tuples.

        Returns:
            None

        """
        self.cds = IntervalTree.from_tuples(cds)

    def set_tid(self, tid: str) -> None:
        """
        Set the transcript ID of the Transcript.

        Args:
            tid (str): The transcript ID.

        Returns:
            None

        """
        self.tid = tid

    def set_gid(self, gid: str) -> None:
        """
        Set the gene ID of the Transcript.

        Args:
            gid (str): The gene ID.

        Returns:
            None

        """
        self.gid = gid

    def nume(self):
        """
        Get the number of exons in the Transcript.

        Returns:
            int: The number of exons.

        """
        return len(self.exons)

    def numc(self):
        """
        Get the number of CDS regions in the Transcript.

        Returns:
            int: The number of CDS regions.

        """
        return len(self.cds)

    def has_cds(self):
        """
        Check if the Transcript has CDS regions.

        Returns:
            bool: True if the Transcript has CDS regions, False otherwise.

        """
        return len(self.cds) > 0

    def introns_it(self):
        """Yield all introns of the transcript.

        Intervals are also inclusive of the lower bound but non-inclusive of the upper bound.
        [first position of the intron, last position of the intron)

        Yields:
            Interval: An Interval object representing an intron.

        """
        if len(self.exons) > 1:
            prev_exon = None
            for e in sorted(self.exons):
                if prev_exon is not None:
                    yield Interval(prev_exon[1], e[0])
                prev_exon = e

    def get_exons(self):
        """Get the sorted list of exons in the Transcript.

        Returns:
            list: A list of exon intervals represented as tuples.

        """
        return sorted(self.exons)

    def get_cds(self):
        """Get the sorted list of CDS regions in the Transcript.

        Returns:
            list: A list of CDS intervals represented as tuples.

        """
        return sorted(self.cds)

    def get_tid(self):
        """Get the transcript ID of the Transcript.

        Returns:
            str: The transcript ID.

        """
        return self.tid

    def get_gid(self):
        """Get the gene ID of the Transcript.

        Returns:
            str: The gene ID.

        """
        return self.gid
    
    def get_cstart(self):
        """Get the start position of the first CDS region in the Transcript.

        Returns:
            int: The start position of the first CDS region.

        """
        return sorted(self.cds)[0][0]

    def get_cend(self):
        """Get the end position of the last CDS region in the Transcript.

        Returns:
            int: The end position of the last CDS region.

        """
        return sorted(self.cds)[-1][1]

    def to_gtf(self):
        """Convert the Transcript to GTF format.

        Returns:
            str: The Transcript object represented in GTF format.

        """
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
            res += e[2].to_gtf() + "\n"
            
        for c in sorted(self.cds):
            res += c[2].to_gtf() + "\n"

        return res.rstrip("\n")

    def to_gff(self):
        """Convert the Transcript to GFF format.

        Returns:
            str: The Transcript object represented in GFF format.

        """
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
            res += e[2].to_gff() + "\n"
            
        for c in sorted(self.cds):
            res += c[2].to_gff() + "\n"

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
    """
    A class representing an exon.

    Attributes:
        obj_type (Types): The type of the object (Exon).

    Inherits from:
        Object

    """
    def __init__(self):
        super().__init__()
        self.tid = None
        self.gid = None

        self.obj_type = Types.Exon

    def __init__(self,obj: Object):
        """Initialize an Exon object.

        Args:
            obj (Object, optional): An Object from which to initialize the Exon. Defaults to None.

        Raises:
            Exception: If a wrong object type is passed to the constructor.

        """
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
        """Initialize the Exon from an Object.

        Args:
            obj (Object): The Object to initialize the Exon from.

        Raises:
            AssertionError: If a wrong object type is passed to the method.

        """
        assert obj.get_type() == Types.Exon or obj.get_type() ==  Types.CDS,"wrong object type passed to Transcript constructor"
        
        self.attrs = rename_attributes(obj.attrs,{"Parent":"transcript_id"})

        assert "transcript_id" in self.attrs,"transcript_id not found in attributes"
        self.tid = self.attrs.get("transcript_id",None)
        if "gene_id" in self.attrs:
            self.gid = self.attrs.get("gene_id",None)

        self.start = obj.start
        self.end = obj.end

    def set_tid(self, tid: str) -> None:
        """Set the transcript ID.

        Args:
            tid (str): The transcript ID.

        """
        self.tid = tid

    def set_gid(self, gid: str) -> None:
        """Set the gene ID.

        Args:
            gid (str): The gene ID.

        """
        self.gid = gid

    def get_tid(self) -> str:
        """Get the transcript ID.

        Returns:
            str: The transcript ID.

        """
        return self.tid

    def get_gid(self) -> str:
        """Get the gene ID.

        Returns:
            str: The gene ID.

        """
        return self.gid

    def get_cstart(self):
        """Get the start position of the first CDS region in the Exon.

        Returns:
            int: The start position of the first CDS region.

        """
        return sorted(self.cds)[0][0]

    def get_cend(self):
        """Get the end position of the last CDS region in the Exon.

        Returns:
            int: The end position of the last CDS region.

        """
        return sorted(self.cds)[-1][1]

    def to_gtf(self):
        """Convert the Exon to GTF format.

        Returns:
            str: The Exon object represented in GTF format without new line at the end.

        """
        res = self.seqid+"\t"+\
                self.source+"\t"+\
                Types.type2str(self.obj_type) +"\t"+\
                str(self.start)+"\t"+\
                str(self.end-1)+"\t"+\
                "."+"\t"+\
                self.strand+"\t"+\
                "."+"\t"+\
                to_attribute_string(self.attrs,False,"exon")
        return res
    
    def to_gff(self):
        """Convert the Exon to GFF format.

        Returns:
            str: The Exon object represented in GFF format  without new line at the end.

        """
        res = self.seqid+"\t"+\
                self.source+"\t"+\
                Types.type2str(self.obj_type) +"\t"+\
                str(self.start)+"\t"+\
                str(self.end-1)+"\t"+\
                "."+"\t"+\
                self.strand+"\t"+\
                "."+"\t"+\
                to_attribute_string(self.attrs,True,"exon")
        return res
    
    __repr__ = to_gtf

    __str__ = __repr__

class CDS(Exon):
    def __init__(self):
        Exon.__init__(self)
        self.tid = None
        self.gid = None

        self.obj_type = Types.CDS

    def __init__(self,obj: Object):
        """Initialize a CDS object.

        Args:
            obj (Object, optional): An Object from which to initialize the CDS. Defaults to None.

        Raises:
            Exception: If a wrong object type is passed to the constructor.

        """
        Exon.__init__(self)
        self.tid = None
        self.gid = None

        self.obj_type = Types.CDS
        
        if obj.get_type() == Types.CDS:
            self.from_exon(obj)
        else:
            raise Exception("wrong object type passed to CDS constructor")

class GTFObjectFactory:
    """
    Object Factory Class.
    Automatically generates objects of correct type from GTF/GFF lines

    """
    @staticmethod
    def create(line):
        """Create an object based on the GTF line.

        Args:
            line (str): A line from the GTF file.

        Returns:
            Object: The created object based on the GTF line.

        """
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
        


# every GTF object has some shared properties
# every object has coordinate range seqid, strand, start and end.
# every object can not span multiple chromosomes and seqids
# every object has source (otherwise set as default)
# every object can (but does not have to) have attributes