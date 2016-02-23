from Bio.Blast.Applications import NcbiblastpCommandline
from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors
from reportlab.lib.colors import red, grey, orange, green, brown, blue, lightblue, purple, yellow

def get_feature(features, id, tags=["protein_id","gene"]):

    """Search list of SeqFeature objects for an identifier under the given tags."""
    for f in features:
        for key in tags:
            #tag may not be present in this feature
            for x in f.qualifiers.get(key, []):
                #if x == id and f.qualifiers.get("gene") != ['parA']:
                if x == id:
                    return f
    raise KeyError(id)

def getCrossLinks(filename):
    A_v_B = [(1,"gg", "gg" )]#structure crosslink to start off (percentmatch, "gene1", "gene2")
    with open(filename) as f:
    # read the file line by line
        for line in f:
            try:
                words = line.split(',')
                A_v_B.append((float(words[2]), words[0] ,words[1]))
            except IndexError:
                print "empty line"
    return A_v_B


C_colors = [yellow]*1+[orange]*1+[brown]*1+[lightblue]*1+[purple]*1+[green]*1+[grey]*1#this creates an array of color for the arrows in the GUI
i = 0 #index of random color to add
geneColor = grey#color of gene. Grey= no name

fileName = "ProuxFigTest25"#name of GUI file to be created
blastfile = "blastfile.csv"#name of the blast file with the results
gd_diagram = GenomeDiagram.Diagram(fileName)
max_len = 0
A_rec = SeqIO.read("pRSB105.gb", "gb")#1st Genbank file to be added to the GUI
B_rec = SeqIO.read("Rms149.gb", "gb")#2nd Genbank file to be added to the GUI
Gname='nn'#name of gene to add


#First section gets the crosslinks from the blast files

A_vs_B = getCrossLinks(blastfile)

print ('(percent, Gene Query, Gene result)')#This prints the list of Blast results for reference
for item in  A_vs_B:
    print item

#this loop adds each gene feature to the record with a color and name
for record, gene_colors in zip([A_rec, B_rec], [C_colors, C_colors]):

    max_len = max(max_len, len(record))
    gd_track_for_features = gd_diagram.new_track(1,
                            name=record.name,
                            greytrack=True,
                            start=0, end=len(record))
    gd_feature_set = gd_track_for_features.new_set()

    for feature in record.features:
        if feature.type != "gene":
            #Exclude this feature
            continue
        try:
            Gname=feature.qualifiers['gene'][0]
            geneColor=gene_colors[i%6]
        except KeyError:#if no gene name make it grey
            Gname ='No Name'
            geneColor = grey;
        gd_feature_set.add_feature(feature, sigil="ARROW",#this adds gene features to gd_feature_set
                                   arrowhead_length = .25,
                                   color = geneColor, label = True,
                                   name = Gname,
                                   label_position = "start",
                                   label_size = 6, label_angle = 0)
        i+=1#increment i so that arrows will have a random color

track_X = gd_diagram.tracks[2]
track_Y = gd_diagram.tracks[1]

#this loop adds the cross links so they point to their feature in the diagram
for score, id_X, id_Y in A_vs_B:
    try:
        feature_X = get_feature(A_rec.features, id_X)
        feature_Y = get_feature(B_rec.features, id_Y)
        color = colors.linearlyInterpolatedColor(colors.white, colors.firebrick, 0, 100, score)
        link_xy = CrossLink((track_X, feature_X.location.start, feature_X.location.end),
                         (track_Y, feature_Y.location.start, feature_Y.location.end),
                         color, colors.lightgrey)

        gd_diagram.cross_track_links.append(link_xy)
    except KeyError:
        print "Feature qualifier for crosslink not found"# for those pesky nameless genes

gd_diagram.draw(format="linear", pagesize=(1200,2400), fragments=1,start=0, end=max_len)
print max_len
gd_diagram.write(fileName + ".pdf", "PDF")
