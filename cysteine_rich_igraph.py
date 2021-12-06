''' Use cyteine-rich blast results to plot 
    the relateness of the genes in igraph
'''

import re
import math
from optparse import OptionParser
from itertools import combinations
from igraph import *

from path import *
from sources import *

parser = OptionParser("usage: %prog [options] blast.tab")
parser.add_option("-s", "--species", dest="species", default = "spiro",
                  help="species the blast results are from")
parser.add_option("-e", "--evalue", dest="evalue", default = 10, type="float",
                  help="evalue cutoff for display")

(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("Please provide converted blast edges tab output file")

edge_file = args[0]
species = options.species
e_cutoff = options.evalue

spiro_genes = \
    session.query(Gene).\
    filter(Gene.partial == None).\
    filter(Gene.pseudo == False).\
    filter(Gene.description.like("Cysteine-rich %")).\
    all()

giardia_genes = \
    session.query(Eupath).\
    filter(Eupath.organism.like("%wb%")).\
    filter(Eupath.version == "2.5").\
    filter(or_(Eupath.description.like("VSP%"), 
               Eupath.description.like("High cysteine%"))).\
    all()

if species == "spiro":
    genes = spiro_genes
elif species == "giardia":
    genes = giardia_genes
elif species == "both":
    genes = spiro_genes + giardia_genes


geneids = []
gene_class = []
gene_shape = []
gene_size = []
for g in genes:
    geneids.append(g.geneid)
    
    # shape
    if re.match("SS", g.geneid):
        gene_shape.append("circle")
    else:
        gene_shape.append("triangle")
    
    # class    
    if re.match("VSP|Cysteine-rich membrane protein 1", g.description):
        gene_class.append(1)
    elif re.match("High cysteine membrane protein|Cysteine-rich membrane protein 2", g.description):
        gene_class.append(2)
    else:
        gene_class.append(3)
    
    # size
    gene_size.append(math.log(g.size, 8)*1.3)


# Read blast info to a dictionary 
edge_info = defaultdict(dict)
with open (edge_file, 'r') as inh:
    for line in inh:
        if re.match("qry", line) or len(line)==0:
            continue
        
        qry_geneid, hit_geneid, score, evalue, p_identity = line.split()
        score = float(score)
        evalue = float(evalue)
        p_identity = float(p_identity)
        edge_info[qry_geneid][hit_geneid] = [score, evalue, p_identity]            


def evalues_ok(edge_info, id1, id2):
    if id2 in edge_info[id1] and id1 in edge_info[id2]:
        evalue1 = edge_info[id1][id2][1]
        evalue2 = edge_info[id2][id1][1]
        return evalue1 < e_cutoff and evalue2 < e_cutoff
    else:
        return False

def get_score(edge_info, id1, id2):
    ''' old: w = 1 - 1/2*(S(a,b)/S(a,a) + S(b,a)/S(b,b))
        w =  S(a,b)/S(a,a) + S(b,a)/S(b,b) ; 0 < w < 2
    '''
    s_aa = edge_info[id1][id1][0]
    s_ab = edge_info[id1][id2][0] if id2 in edge_info[id1] else 0
    s_ba = edge_info[id2][id1][0] if id1 in edge_info[id2] else 0
    s_bb = edge_info[id2][id2][0]
    return s_ab/float(s_aa) + s_ba/float(s_bb)

def edge_weight(score):
    return (10 ** score - 1) / 90 / 3.3 # score / 8
    

def edge_color(score):
    return "#000000%02x" % int((score * 127.5/3.3))

def add_edge(g, i, j, score):
    g.add_edges([(i,j)])
    edge_id = g.get_eid(i,j)
    weight = edge_weight(score)
    g.es[edge_id]["weight"] = weight
    g.es[edge_id]['color'] = edge_color(score)


g = Graph()
g.add_vertices(len(geneids))
g.vs["name"] = geneids
g.vs["class"] = gene_class
g.vs['shape'] = gene_shape
color_dict = {1: "blue", 2: "pink", 3: "green"}
g.vs["color"] = [color_dict[cl] for cl in g.vs["class"]]
g.vs["size"] = gene_size

scores = []
for i, j in combinations(range(len(geneids)), 2):
    id1, id2 = geneids[i], geneids[j]
    if evalues_ok(edge_info, id1, id2):
        score = get_score(edge_info, id1, id2)
        add_edge(g, i, j, score)
        scores.append(score)

# outfile = "plots/igraph/%s_crp.w2.%.0e.%s.pdf" %(species, e_cutoff, "drl")
# layout = g.layout_drl(weights = g.es["weight"])
# plot(g, target = outfile, layout = layout, bbox = (600, 600), margin = 50, vertex_size=7)

outfile = "plots/igraph/%s_crp.w2.%.0e.%s.pdf" %(species, e_cutoff, "fr")
layout = g.layout_fruchterman_reingold(weights=g.es["weight"])
# plot(g, target =outfile, layout = layout, bbox = (600, 600), margin = 30)
plot(g, target =outfile, layout = layout, bbox = (260, 260), margin = 5) # 3.6 inches


from rpy2.robjects import *

r('weights <- %s' % FloatVector(g.es["weight"]).r_repr()) 
r('scores <- %s' % FloatVector(scores).r_repr()) 
r('''
pdf(file = "plots/igraph/igraph_weights.pdf", width = 3.15, height = 1.75, pointsize = 8)
hist(scores)
hist(weights)
dev.off()
''')
                