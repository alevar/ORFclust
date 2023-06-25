import unittest

from classes.transcript import GTFObjectFactory, Object, Transcript, Exon, CDS

class MyTestCase(unittest.TestCase):
    def test_create_transcript_1(self):
        line = """chr1    ORFclust        transcript      17308197        17356057        .       +       .       transcript_id "rna-XM_011541154.3"; gene_id "gene-PADI4"; gene_name "PADI4"; description "peptidyl arginine deiminase 4"; CDS_Dbxref "GeneID:23569,Genbank:XP_011539456.1,HGNC:HGNC:18368,MIM:605347"; CDS_Name "XP_011539456.1"; CDS_gbkey "CDS"; CDS_product "protein-arginine deiminase type-4 isoform X4"; Dbxref "GeneID:23569,Genbank:XM_011541154.3,HGNC:HGNC:18368,MIM:605347"; Name "XM_011541154.3"; gbkey "mRNA"; gene "PADI4"; gene_biotype "protein_coding"; gene_synonym "PAD,PAD4,PADI5,PDI4,PDI5"; model_evidence "Supporting evidence includes similarity to: 17 ESTs%2C 1 long SRA read%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 3 samples with support for all annotated introns"; product "peptidyl arginine deiminase 4%2C transcript variant X5"; protein_id "XP_011539456.1"; transcriptID "XM_011541154.3";"""
        
        obj = Object()
        obj.add_line(line)
        tx = Transcript(obj)
        self.assertEqual(True, True)

if __name__ == '__main__':
    unittest.main()
