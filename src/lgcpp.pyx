import sys
from math import log
from collections import defaultdict
from cython.operator cimport dereference as deref, preincrement as inc
from cython import address as addr

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.map cimport map

cimport cython

ctypedef fused int_or_str:
    cython.int
    cython.p_char
    bytes
    unicode
    string
    
cdef extern from "math.h": 
    bint isnan(double x)

cdef extern from "superdouble.h":
    cdef cppclass _Superdouble "Superdouble":
        _Superdouble()
        _Superdouble(long double, int)
        _Superdouble operator/ (_Superdouble)
        _Superdouble operator+ (_Superdouble)
        _Superdouble operator- (_Superdouble)
        void operator++ ()
        void operator-- ()
        ## void operator*= (Superdouble)
        ## void operator/= (Superdouble)
        ## void operator+= (Superdouble)
        ## void operator-= (Superdouble)
        bool operator < (_Superdouble&)
        bool operator > (_Superdouble&)
        bool operator >= (_Superdouble&)
        bool operator <= (_Superdouble&)
        bool operator == (_Superdouble&)
        int getExponent()
        _Superdouble add(_Superdouble)
        void divideby(_Superdouble)
        long double getMantissa()
        _Superdouble getLn()
        _Superdouble abs()
        void switch_sign()
        void adjustDecimal()

    _Superdouble super_add(_Superdouble, _Superdouble)
    _Superdouble super_divide(_Superdouble, _Superdouble)
    _Superdouble super_ln(_Superdouble)

cdef _Superdouble * superptr(_Superdouble x):
    return &x

cdef class Superdouble

cdef class Superdouble:
    cdef _Superdouble ptr
    def __cinit__(self, long double mantissa=1.0, int exponent=0):
        self.ptr = _Superdouble(mantissa, exponent)
    @property
    def mantissa(self):
        return self.ptr.getMantissa()
    @property
    def exponent(self):
        return self.ptr.getExponent()
    def __str__(self):
        self.ptr.adjustDecimal()
        return "%ge%d" % (self.ptr.getMantissa(), self.ptr.getExponent())
    def __repr__(self):
        self.ptr.adjustDecimal()
        return "Superdouble(mantissa={}, exponent={})".format(
            self.ptr.getMantissa(), self.ptr.getExponent())
    def __float__(self):
        return float(str(self))
    def __add__(Superdouble self, Superdouble other):
        return superdouble_factory(super_add(self.ptr, other.ptr))
    def __truediv__(Superdouble self, Superdouble other):
        cdef long double m = self.ptr.getMantissa()/other.ptr.getMantissa()
        cdef int e = self.ptr.getExponent()-other.ptr.getExponent()
        return Superdouble(m, e)
    def __eq__(self, Superdouble other):
        return self.ptr == other.ptr
    def __gt__(self, Superdouble other):
        return self.ptr > other.ptr
    def __ge__(self, Superdouble other):
        return self.ptr >= other.ptr
    def __lt__(self, Superdouble other):
        return self.ptr < other.ptr
    def __le__(self, Superdouble other):
        return self.ptr <= other.ptr
    def getLn(self):
        return float(superdouble_factory(super_ln(self.ptr)))
    def adjustDecimal(self):
        self.ptr.adjustDecimal()


cdef Superdouble superdouble_factory(_Superdouble p):
    cdef Superdouble n = Superdouble.__new__(Superdouble)
    p.adjustDecimal()
    n.ptr = p
    return n

## cdef long double super2double(_Superdouble x):
##     return x.getMantissa()*pow(10., x.getExponent())

cdef extern from "BranchSegment.h":
    cdef cppclass _BranchSegment "BranchSegment":
        _BranchSegment(double, int)
        double duration
        int period
        int startdistint
        _RateModel * model
        vector[_Superdouble]* distconds
        vector[int] fossilareaindices
        vector[int] getFossilAreas()
        void setFossilArea(int)
        void setModel(_RateModel*)
        
cdef class BranchSegment:
    cdef _BranchSegment* ptr
    def __init__(self, double duration, int period, RateModel m):
        cdef _BranchSegment * bs = new _BranchSegment(duration, period)
        self.ptr = bs
        self.ptr.setModel(m.ptr)
    def __cinit__(self):
        self.ptr = NULL
    ## def __dealloc__(self):
    ##     del self.ptr
    property duration:
        def __get__(self):
            return self.ptr.duration
        def __set__(self, double x):
            self.ptr.duration = x
    property period:
        def __get__(self):
            return self.ptr.period
        def __set__(self, int x):
            self.ptr.period = x
    ## property fossilareaindices:
    ##     def __get__(self):
    ##         return self.ptr.fossilareaindices
    def setFossilArea(self, int area):
        ## print >> sys.stderr, 'setFossilArea: %s' % area
        self.ptr.setFossilArea(area)
    def getFossilAreas(self):
        return self.ptr.getFossilAreas()
    def setModel(self, RateModel m):
        self.ptr.setModel(m.ptr)

cdef BranchSegment branchsegment_factory(_BranchSegment *p):
    cdef BranchSegment bs = BranchSegment.__new__(BranchSegment)
    bs.ptr = p
    return bs

def new_branchsegment(double duration, int period, RateModel m):
    cdef _BranchSegment * newbs = new _BranchSegment(duration, period)
    cdef BranchSegment bs = BranchSegment.__new__(BranchSegment)
    bs.ptr = newbs
    bs.setModel(m)
    return bs

cdef extern from "node.h":
    cdef cppclass _Node "Node":
        _Node()

        double BL
        double height
        int number
        string name
        _Node * parent
        vector[_Node*] children
        string comment
        vector[_BranchSegment]* segs
        vector[vector[int]]* excluded_dists

        bool isExternal()
        bool isInternal()
        bool isRoot()
        vector[_Node*] getChildren()
        string getName()
        vector[_BranchSegment]* getSegVector()
        string getNewick(bool)
        string getNewick(bool,string)
        int getNumber()
        _Node * getParent()
        vector[vector[int]] * getExclDistVector()
        double getHeight()
        double getBL()
        double lengthToRoot()
        void setHeight(double h)
        void deleteSegVector()
        void initSegVector()
        
cdef class Node
cdef class Node:
    cdef _Node* ptr
    def __cinit__(self):
        self.ptr = NULL
    ## def __dealloc__(self):
    ##     del self.ptr
    def isExternal(self):
        return self.ptr.isExternal()
    def isInternal(self):
        return self.ptr.isInternal()
    def isRoot(self):
        return self.ptr.isRoot()
    def getChildren(self):
        cdef vector[_Node*] v = self.ptr.getChildren()
        cdef vector[_Node*].iterator it = v.begin()
        children = []
        while it != v.end():
            children.append(node_factory(deref(it)))
            inc(it)
        return children
    def iterChildren(self):
        cdef vector[_Node*] v = self.ptr.getChildren()
        cdef vector[_Node*].iterator it = v.begin()
        while it != v.end():
            yield node_factory(deref(it))
            inc(it)
    def getName(self):
        return self.ptr.getName().decode('utf-8')
    def getBL(self):
        return self.ptr.getBL()
    def getNumber(self):
        cdef int n = self.ptr.getNumber()
        return n
    def getHeight(self):
        return self.ptr.getHeight()
    def setHeight(self, double height):
        self.ptr.setHeight(height)
    def getParent(self):
        cdef _Node* p = self.ptr.getParent()
        if p != NULL:
            return node_factory(self.ptr.getParent())
        return None
    def getSegVector(self):
        cdef vector[_BranchSegment] *v = self.ptr.getSegVector()
        cdef vector[_BranchSegment].iterator it = deref(v).begin()
        segs = []
        while it != deref(v).end():
            segs.append(branchsegment_factory(&deref(it)))
            inc(it)
        return segs
    def lengthToRoot(self):
        return self.ptr.lengthToRoot()

    ## def set_tip_conditional(self, double i):
    ##     cdef vector[_BranchSegment]* v = self.ptr.getSegVector()
    ##     cdef _BranchSegment seg = deref(v.at(0))
    ##     seg.distconds.at(i) = <double>i

    def getExclDistVector(self):
        return deref(self.ptr.getExclDistVector())

    def iternodes(self, f=None):
        """
        generate nodes descendant from self - including self
        """
        if (f and f(self)) or (not f):
            yield self
        for child in self.iterChildren():
            for n in child.iternodes():
                if (f and f(n)) or (not f):
                    yield n

    def preiter(self, f=None):
        for n in self.iternodes(f=f):
            yield n

    def rootpath_length(self):
        cdef double bl = 0
        cdef Node n = self
        while not n.ptr.isRoot():
            bl += n.ptr.getBL()
            n = n.getParent()
        return bl

    # def setSegVector(self, durations, periods, models):
    #     self.ptr.deleteSegVector()
    #     cdef vector[_BranchSegment]* segs = new vector[_BranchSegment]()
    #     cdef BranchSegment bs
    #     v = []
    #     for d, p, m in zip(durations, periods, models):
    #         bs = BranchSegment(d, p)
    #         bs.setModel(m)
    #         segs.push_back(deref(bs.ptr))
    #         v.append(bs)
    #     self.ptr.segs = segs
    #     return v

    def setSegVector(self, v):
        self.ptr.deleteSegVector()
        self.ptr.initSegVector()
        ## cdef vector[_BranchSegment]* segs = new vector[_BranchSegment]()
        cdef BranchSegment bs
        for bs in v:
            self.ptr.segs.push_back(deref(bs.ptr))
        ## self.ptr.segs = segs


cdef Node node_factory(_Node *p):
    cdef Node n = Node.__new__(Node)
    n.ptr = p
    return n

cdef extern from "RateModel.h":
    cdef cppclass _RateModel "RateModel":
        _RateModel(int, bool, vector[double], bool)
        vector[string] get_labels()
        void set_nthreads(int)
        int get_nthreads()
        void setup_dists()
        void setup_dists(vector[vector[int]], bool)
        void setup_Dmask()
        void set_Dmask_cell(int, int, int, double, bool)
        void setup_D(double)
        void setup_D_provided(double, vector[vector[vector[double]]])
        void setup_E(double)
        void set_Qdiag(int)
        void setup_Q()
        vector[vector[int]] * getDists()
        
cdef class RateModel:
    cdef _RateModel* ptr
    def __cinit__(self, int na, bool ge, list periods, bool is_sparse):
        cdef vector[double] v = vector[double]()
        for x in periods: v.push_back(x)
        self.ptr = new _RateModel(na, ge, v, is_sparse)
    ## def __dealloc__(self):
    ##     del self.ptr
    def set_nthreads(self, int n):
        self.ptr.set_nthreads(n)
    def get_nthreads(self):
        return self.ptr.get_nthreads()
    def setup_dists(self, indists=[], bool incl=True):
        cdef vector[vector[int]] v = vector[vector[int]]()
        cdef vector[int]* k
        if len(indists):
            for row in indists:
                k = new vector[int]()
                for x in row: k.push_back(x)
                v.push_back(deref(k))
            self.ptr.setup_dists(v, incl)
        else:
            self.ptr.setup_dists()
    def setup_D(self, double x):
        self.ptr.setup_D(x)
    def setup_D_provided(self, double d, vector[vector[vector[double]]] dmask):
        self.ptr.setup_D_provided(d, dmask)
    def setup_E(self, double x):
        self.ptr.setup_E(x)
    def setup_Q(self):
        self.ptr.setup_Q()
    def setup_Dmask(self):
        self.ptr.setup_Dmask()
    def set_Dmask_cell(self, int period, int area, int area2,
                       double prob, bool sym):
        self.ptr.set_Dmask_cell(period, area, area2, prob, sym)
    def getDists(self):
        cdef vector[vector[int]] * dists = self.ptr.getDists()
        return deref(dists)

cdef extern from "tree.h":
    cdef cppclass _Tree "Tree":
        _Tree()
        int getNodeCount()
        int getExternalNodeCount()
        int getInternalNodeCount()
        _Node* getExternalNode(string)
        _Node* getExternalNode(int)
        _Node* getInternalNode(int)
        _Node* getInternalNode(string &)
        _Node* getRoot()
        _Node* getMRCA(vector[string] innodes)
        double maxTipPathLength()

cdef class Tree:
    cdef _Tree* ptr
    def __cinit__(self):
        self.ptr = new _Tree()
    ## def __dealloc__(self):
    ##     del self.ptr
    def getNodeCount(self):
        return self.ptr.getNodeCount()
    def getExternalNodeCount(self):
        return self.ptr.getExternalNodeCount()
    def getInternalNodeCount(self):
        return self.ptr.getInternalNodeCount()
    def getExternalNode(self, string s):
        cdef _Node* p = self.ptr.getExternalNode(s)
        return node_factory(p)
    ## def getExternalNode(self, int i):
    ##     cdef _Node* p = self.ptr.getExternalNode(i)
    ##     return node_factory(p)
    def getInternalNode(self, int_or_str x):
        cdef _Node* p
        cdef Node n
        cdef string s
        if int_or_str is cython.int:
            ## print >> sys.stderr, 'getInternalNode got int', x
            p = self.ptr.getInternalNode(x)
        elif int_or_str is unicode:
            ## print >> sys.stderr, 'getInternalNode got str', x
            s = x.encode('utf-8')
            p = self.ptr.getInternalNode(s)
        else:
            raise ValueError
        if p != NULL:
            return node_factory(p)
        else:
            raise Exception
    def getRoot(self):
        cdef _Node* p = self.ptr.getRoot()
        return node_factory(p)
    def getMRCA(self, names):
        cdef vector[string] v = vector[string]()
        for s in names: v.push_back(string(<char *>s))
        cdef _Node* p = self.ptr.getMRCA(v)
        return node_factory(p)
    def internalNodes(self):
        cdef int n = self.ptr.getInternalNodeCount()
        return [ node_factory(self.ptr.getInternalNode(i)) for i in range(n) ]
    def setHeights(self):
        cdef double maxh
        cdef Node root = self.getRoot()
        cdef Node n
        cdef _Node* nptr
        leaves = root.iternodes(lambda x:x.isExternal())
        maxh = max([ lf.rootpath_length() for lf in leaves ])
        root.ptr.setHeight(maxh)
        it = root.iternodes()
        next(it)
        for n in it:
            nptr = n.ptr
            nptr.setHeight(nptr.getParent().getHeight()-nptr.getBL())
    def maxTipPathLength(self):
        return self.ptr.maxTipPathLength()
    def newick(self):
        #cdef string s = string(<char *>"number")
        return "".join([self.ptr.getRoot().getNewick(True).decode('utf-8'),';'])

cdef extern from "tree_reader.h":
    cdef cppclass _TreeReader "TreeReader":
        _TreeReader()
        _Tree* readTree(string)

def readtree(s):
    cdef _TreeReader* reader = new _TreeReader()
    cdef string treestr = s.encode('utf-8')
    cdef _Tree* tree = reader.readTree(<string>treestr)
    cdef Tree t = Tree()
    t.ptr = tree
    ## del reader
    return t

cdef extern from "AncSplit.h":
    cdef cppclass _AncSplit "AncSplit":
        _AncSplit(_RateModel*, int, int, int, _Superdouble)
        _RateModel* getModel()
        double getWeight()
        _Superdouble getLikelihood()
        int ancdistint, ldescdistint, rdescdistint

## cdef class AncSplit:
##     cdef _AncSplit* ptr
##     def __cinit__(self, RateModel m, int dist, int ldesc, int rdesc, Superdouble w):
##         self.ptr = new _AncSplit(m.ptr, dist, ldesc, rdesc, deref(w.ptr))

cdef extern from "OptimizeBioGeo.h":
    cdef cppclass _OptimizeBioGeo "OptimizeBioGeo":
        _OptimizeBioGeo(_BioGeoTree*, _RateModel*, bool)
        vector[double] optimize_global_dispersal_extinction()

cdef extern from "BioGeoTree.h":
    cdef cppclass _BioGeoTree "BioGeoTree":
        _BioGeoTree(_Tree *, vector[double])
        void set_default_model(_RateModel *)
        _Superdouble eval_likelihood(bool)
        void set_tip_conditionals(map[string,vector[int]])
        void ancdist_conditional_lh(_Node &, bool)
        void set_store_p_matrices(bool)
        void set_use_stored_matrices(bool)
        void update_default_model(_RateModel * mod)
        void prepare_ancstate_reverse()
        map[vector[int],vector[_AncSplit]] calculate_ancsplit_reverse(_Node &,bool)
        vector[_Superdouble] calculate_ancstate_reverse(_Node &,bool)
        void set_excluded_dist(vector[int], _Node*)
        void setFossilatNodeByMRCA(vector[string], int)
        void setFossilatNodeByMRCA_id(_Node *, int)
        void setFossilatBranchByMRCA(vector[string], int, double)
        void setFossilatBranchByMRCA_id(_Node *, int, double)
        
cdef class BioGeoTree:
    cdef _BioGeoTree* ptr
    def __cinit__(self, Tree t, list periods):
        cdef vector[double] v = vector[double]()
        for x in periods: v.push_back(x)
        self.ptr = new _BioGeoTree(t.ptr, v)
    ## def __dealloc__(self):
    ##     del self.ptr
    def set_store_p_matrices(self, bool b):
        self.ptr.set_store_p_matrices(b)
    def set_use_stored_matrices(self, bool b):
        self.ptr.set_use_stored_matrices(b)
    def update_default_model(self, RateModel m):
        self.ptr.update_default_model(m.ptr)
    def set_default_model(self, RateModel m):
        self.ptr.set_default_model(m.ptr)
    def set_excluded_dist(self, dist, Node n):
        cdef vector[int] v = vector[int]()
        cdef int x
        for x in dist: v.push_back(x)
        ## print >> sys.stderr, 'excluding {} from node {}'.format(
        ##     dist, n.getName())
        self.ptr.set_excluded_dist(v, n.ptr)
    def set_tip_conditionals(self, data):
        cdef map[string,vector[int]] m #= map[string,vector[int]]()
        #cdef string* s
        cdef vector[int]* dist
        cdef string tiplabel
        for k, v in sorted(data.items()):
            #s = new string(<char *>k)
            tiplabel = k.encode('utf-8')
            dist = new vector[int]()
            for x in v:
                dist.push_back(x)
            m[tiplabel] = deref(dist)
        self.ptr.set_tip_conditionals(m)

    def optimize_global_dispersal_extinction(self, bool marginal, RateModel m):
        ## cdef double initL = super2double(self.ptr.eval_likelihood(marginal))
        print >> sys.stderr, "optimizing rate parameters..."
        cdef _OptimizeBioGeo* opt = new _OptimizeBioGeo(self.ptr, m.ptr, marginal)
        cdef vector[double] disext = opt.optimize_global_dispersal_extinction()
        cdef double d, e, neglnL
        d = disext[0]; e = disext[1]
        ## print "dispersal rate:", disext[0]
        ## print "local extinction rate:", disext[1]
        ## print
        m.setup_D(disext[0])
        m.setup_E(disext[1])
        m.setup_Q()
        self.ptr.update_default_model(m.ptr)
        self.ptr.set_store_p_matrices(True)
        neglnL = float(superdouble_factory(self.ptr.eval_likelihood(marginal)))
        ## neglnL = float(superdouble_factory(superptr(self.ptr.eval_likelihood(marginal))))
        ## print "-lnL:", neglnL
        ## print
        ## self.ptr.set_store_p_matrices(False)
        print >> sys.stderr, ("dispersal = %s; extinction = %s; -lnL = %s" %
                              (d, e, neglnL))
        return (d, e, neglnL)
        
    def set_global_dispersal_extinction(self, double d, double e, bool marginal, RateModel m):
        ## cdef double initL = super2double(self.ptr.eval_likelihood(marginal))
        cdef neglnL
        m.setup_D(d)
        m.setup_E(e)
        m.setup_Q()
        self.ptr.update_default_model(m.ptr)
        self.ptr.set_store_p_matrices(True)
        neglnL = float(superdouble_factory(self.ptr.eval_likelihood(marginal)))
        ## neglnL = float(superdouble_factory(superptr(self.ptr.eval_likelihood(marginal))))
        print >> sys.stderr, ("dispersal = %s; extinction = %s; -lnL = %s" %
                              (d, e, neglnL))

    def ancsplits(self, Tree intree, bool marginal, RateModel m, list areas):
        print >> sys.stderr, "calculating ancestral splits..."
        cdef int n = intree.ptr.getInternalNodeCount()
        #print >> sys.stderr, "%s nodes" % n
        cdef int i, j, k
        cdef _Node* node
        cdef map[vector[int],vector[_AncSplit]] ras
        cdef map[vector[int],vector[_AncSplit]].iterator rasit
        cdef vector[_AncSplit]* tans
        cdef vector[_AncSplit].iterator tansit
        cdef _AncSplit* ancsplit
        cdef map[int,string] area_i2s = map[int,string]()
        cdef _BioGeoTreeTools tt
        cdef map[_Superdouble,string] summary
        cdef map[_Superdouble,string].iterator it
        cdef _Superdouble totL, zero, tmp, prop
        cdef Superdouble lh

        self.ptr.set_use_stored_matrices(True)
        self.ptr.prepare_ancstate_reverse()

        for i, a in enumerate(areas):
            area_i2s[i] = a.encode('utf-8')

        dists = [ tuple(x) for x in m.getDists() ]

        def d2s(d, sep=''):
            return sep.join([ areas[i] for i,x in enumerate(d) if x ])

        zero = _Superdouble(0,0)

        d = defaultdict(list)
        for i in range(n):
            node = intree.ptr.getInternalNode(i)
            name = node.getName().decode('utf-8')
            ## _excl = node.getExclDistVector()
            ## excl = set()
            ## for j in range(_excl.size()):
            ##     excl.add(tuple(_excl.at(j)))
            totL = _Superdouble(0,0)
            ras = self.ptr.calculate_ancsplit_reverse(deref(node), marginal)
            rasit = ras.begin()
            while rasit != ras.end():
                tans = &deref(rasit).second
                tansit = tans.begin()
                while tansit != tans.end():
                    ancsplit = &deref(tansit)
                    tmp = ancsplit.getLikelihood()
                    if tmp > zero:
                        totL = super_add(totL, tmp)
                        d[name].append(
                            [superdouble_factory(tmp),
                             d2s(dists[ancsplit.ancdistint]),
                             d2s(dists[ancsplit.ldescdistint]),
                             d2s(dists[ancsplit.rdescdistint])])
                    inc(tansit)
                inc(rasit)
        print >> sys.stderr, "Done"
        return d

    def setFossilatNodebyMRCA(self, names, int area):
        cdef vector[string] v = vector[string]()
        for s in names: v.push_back(string(<char *>s))
        self.ptr.setFossilatNodeByMRCA(v, area)

    def setFossilatNodebyMRCA_id(self, Node n, int area):
        self.ptr.setFossilatNodeByMRCA_id(n.ptr, area)

    def setFossilatBranchbyMRCA(self, names, int area, double age):
        cdef vector[string] v = vector[string]()
        for s in names: v.push_back(string(<char *>s))
        self.ptr.setFossilatBranchByMRCA(v, area, age)

    def setFossilatBranchbyMRCA_id(self, Node n, int area, double age):
        self.ptr.setFossilatBranchByMRCA_id(n.ptr, area, age)

            
cdef extern from "BioGeoTreeTools.h":
    cdef cppclass _BioGeoTreeTools "BioGeoTreeTools":
        map[_Superdouble,string] summarizeSplits(
            _Node *,
            map[vector[int],vector[_AncSplit]] &,
            map[int,string] &, _RateModel *
            )
        
cdef extern from "InputReader.h":
    cdef cppclass _InputReader "InputReader":
        _InputReader()
        void readMultipleTreeFile(string, vector[_Tree*])
        map[string,vector[int]] readStandardInputData(string)
        void checkData(map[string,vector[int]],vector[_Tree])
        int nareas
        int nspecies

cdef class InputReader:
    cdef _InputReader* ptr
    ## cdef vector[_Tree*]* trees
    def __cinit__(self):
        self.ptr = new _InputReader()
        ## self.trees = new vector[_Tree*]()
    ## def __dealloc__(self):
    ##     del self.ptr

    def read_treefile(self, filename):
        cdef string s = string(<char *>filename)
        cdef vector[_Tree*] *trees = new vector[_Tree*]() #deref(self.trees)
        cdef vector[_Tree*] tv
        self.ptr.readMultipleTreeFile(<string>s, deref(trees))
        tv = deref(trees)
        print tv.size()
