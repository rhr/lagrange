import lgcpp, ivy
from pprint import pprint
from ivy.interactive import *

def rng2dist(rng, nareas):
    v = [0]*nareas
    for i in rng: v[i] = 1
    return v

with open("test.data") as f:
    #m, t, d, nodelabels, base_rates = lagrange.input.eval_decmodel(f.read())
    d = eval(f.read())

areas = d['area_labels']
nareas = len(areas)
ranges = d['ranges']
data = d['taxon_range_data']
newick = d['newick_trees'][0]['newick']

def showtree(r):
    from ivy.vis.symbols import tipsquares
    tango = ivy.vis.colors.tango()
    gray = ivy.vis.colors.tango_colors['Aluminium3']
    colors = [ tango.next() for x in areas ]

    f = treefig(r)
    leaves = f.root.leaves()
    for lf in leaves:
        v = [gray]*nareas
        for i, val in enumerate(data[lf.label]):
            if val: v[i] = colors[i]
        f.detail.decorate(tipsquares, lf, v, size=8)

V = [ [ 0 for x in areas ] for r in ranges ]
for i, r in enumerate(ranges):
    for j in r: V[i][j] = 1
    
periods = d['dispersal_durations']

for k, v in data.items():
    data[k] = rng2dist(v, nareas)

t = lgcpp.readtree(newick)
r = ivy.tree.read(t.newick())

ge = True; sparse = False
model = lgcpp.RateModel(len(areas), ge, periods, sparse)
model.set_nthreads(1)
model.setup_Dmask()
model.setup_dists(V)
model.setup_D(0.01)
model.setup_E(0.01)
model.setup_Q()

bgt = lgcpp.BioGeoTree(t, periods)
bgt.set_default_model(model)
bgt.set_tip_conditionals(data)
marginal = True
bgt.optimize_global_dispersal_extinction(marginal, model)
n2split = bgt.ancsplits(t, marginal, model, areas)

#pprint(n2split)

for n in r.preiter(lambda x:x.children):
    i = int(n.label)+1
    print "ancestral splits for node %s:" % i
    for lnl, s in n2split[i]:
        print " ", lnl, s
