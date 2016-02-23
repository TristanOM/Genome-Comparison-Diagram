from Bio.Blast.Applications import NcbiblastpCommandline
from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def writeFasta(Filename, FastaFileName):
    handle1=open(Filename,"rU")
    record1 = SeqIO.read(handle1, "gb")
    outhandle1 = open(FastaFileName,'w')
    for feature in record1.features:
        if feature.type=='CDS':
            proteinseq = feature.extract(record1)
            try:
                proteinseq.id=feature.qualifiers['protein_id'][0]
            except KeyError:#if no protein id to use
                try:
                    print 'protein_id qualifier not found. Using gene qualifier'
                    print feature.qualifiers['gene'][0]
                    proteinseq.id=feature.qualifiers['gene'][0]
                except KeyError:#if no product name
                    proteinseq.id='No_id'
                    print "gene qualifier found. Using No_id as identifier"
            SeqIO.write(proteinseq,outhandle1,"fasta")
    outhandle1.close()

GBFilename2 = "Rms149.gb"
GBFilename1 = "pRSB105.gb"
FastaFileName1 = "file1.fasta"
FastaFileName2 = "file2.fasta"
print('File2')
writeFasta(GBFilename1,FastaFileName1)
print('File1')
writeFasta(GBFilename2,FastaFileName2)