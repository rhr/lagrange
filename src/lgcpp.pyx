from cython.operator cimport dereference as deref, preincrement as inc
from cython import address as addr

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.map cimport map

cdef extern from "superdouble.h":
    cdef cppclass _Superdouble "Superdouble":
        _Superdouble(long double, int)
        _Superdouble operator/ (_Superdouble)
        _Superdouble operator+ (_Superdouble)
        _Superdouble operator- (_Superdouble)
        void operator++ ()
        void operator -- ()
        ## void operator*= (Superdouble)
        ## void operator/= (Superdouble)
        ## void operator+= (Superdouble)
        ## void operator-= (Superdouble)
        ## bool operator < (const Superdouble&)const 
        ## bool operator > (const Superdouble&)const 
        ## bool operator >= (const Superdouble&)const 
        ## bool operator <= (const Superdouble&)const 
        int getExponent()
        double getMantissa()
        _Superdouble getLn()
        _Superdouble abs()
        void switch_sign()

cdef class Superdouble:
    cdef _Superdouble* ptr
    def __cinit__(self, long double mantissa=1.0, int exponent=0):
        self.ptr = new _Superdouble(mantissa, exponent)

cdef double super2double(_Superdouble x):
    return x.getMantissa()*pow(10., x.getExponent())

cdef extern from "node.h":
    cdef cppclass _Node "Node":
        _Node()
        string getName()
        vector[_BranchSegment]* getSegVector()

cdef class Node:
    cdef _Node* ptr
    def __cinit__(self):
        self.ptr = NULL
    def __dealloc__(self):
        del self.ptr
    def getName(self):
        return self.ptr.getName().c_str()
    ## def set_tip_conditional(self, double i):
    ##     cdef vector[_BranchSegment]* v = self.ptr.getSegVector()
    ##     cdef _BranchSegment seg = deref(v.at(0))
    ##     seg.distconds.at(i) = <double>i

cdef Node node_factory(_Node *p):
    cdef Node n = Node.__new__(Node)
    n.ptr = p
    return n

cdef extern from "BranchSegment.h":
    cdef cppclass _BranchSegment "BranchSegment":
        _BranchSegment(double, int)
        vector[_Superdouble]* distconds

cdef class BranchSegment:
    cdef _BranchSegment* ptr
    def __cinit__(self):
        self.ptr = NULL
    def __dealloc__(self):
        del self.ptr
    ## property distconds:
    ##     def __get__(self):
    ##         return self.ptr.distconds

cdef BranchSegment branchsegment_factory(_BranchSegment *p):
    cdef BranchSegment bs = BranchSegment.__new__(BranchSegment)
    bs.ptr = p
    return bs

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
    def __dealloc__(self):
        del self.ptr
    def set_nthreads(self, int n):
        self.ptr.set_nthreads(n)
    def get_nthreads(self):
        return self.ptr.get_nthreads()
    def setup_dists(self, list indists=None, bool incl=True):
        cdef vector[vector[int]] v = vector[vector[int]]()
        cdef vector[int]* k
        if indists:
            for row in indists:
                k = new vector[int]()
                for x in row: k.push_back(x)
                v.push_back(deref(k))
            self.ptr.setup_dists(v, incl)
        else:
            self.ptr.setup_dists()
    def setup_D(self, double x):
        self.ptr.setup_D(x)
    def setup_E(self, double x):
        self.ptr.setup_E(x)
    def setup_Q(self):
        self.ptr.setup_Q()
                
    def setup_Dmask(self):
        self.ptr.setup_Dmask()
    def set_Dmask_cell(self, int period, int area, int area2,
                       double prob, bool sym):
        self.ptr.set_Dmask_cell(period, area, area2, prob, sym)

cdef extern from "tree.h":
    cdef cppclass _Tree "Tree":
        _Tree()
        int getNodeCount()
        int getExternalNodeCount()
        int getInternalNodeCount()
        _Node* getExternalNode(int)
        _Node* getInternalNode(int)

cdef class Tree:
    cdef _Tree* ptr
    def __cinit__(self):
        self.ptr = new _Tree()
    def __dealloc__(self):
        del self.ptr
    def getNodeCount(self):
        return self.ptr.getNodeCount()
    def getExternalNodeCount(self):
        return self.ptr.getExternalNodeCount()
    def getExternalNode(self, int i):
        cdef _Node* p = self.ptr.getExternalNode(i)
        return node_factory(p)
    def getInternalNode(self, int i):
        cdef _Node* p = self.ptr.getInternalNode(i)
        return node_factory(p)
    def internalNodes(self):
        cdef int n = self.ptr.getInternalNodeCount()
        return [ node_factory(self.ptr.getInternalNode(i)) for i in range(n) ]

cdef extern from "tree_reader.h":
    cdef cppclass _TreeReader "TreeReader":
        _TreeReader()
        _Tree* readTree(string)

def readtree(s):
    cdef _TreeReader* reader = new _TreeReader()
    cdef string treestr = string(<char *>s)
    cdef _Tree* tree = reader.readTree(<string>treestr)
    cdef Tree t = Tree()
    t.ptr = tree
    del reader
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

cdef extern from "BioGeoTreeTools.h":
    cdef cppclass _BioGeoTreeTools "BioGeoTreeTools":
        void summarizeSplits(_Node *,
                             map[vector[int],vector[_AncSplit]] &,
                             map[int,string] &, _RateModel *)
        

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
        
cdef class BioGeoTree:
    cdef _BioGeoTree* ptr
    def __cinit__(self, Tree t, list periods):
        cdef vector[double] v = vector[double]()
        for x in periods: v.push_back(x)
        self.ptr = new _BioGeoTree(t.ptr, v)
    def __dealloc__(self):
        del self.ptr
    def set_default_model(self, RateModel m):
        self.ptr.set_default_model(m.ptr)
    def set_tip_conditionals(self, data):
        cdef map[string,vector[int]] m = map[string,vector[int]]()
        #cdef string* s
        cdef vector[int]* dist
        for k, v in sorted(data.items()):
            #s = new string(<char *>k)
            dist = new vector[int]()
            for x in v: dist.push_back(x)
            m[string(<char *>k)] = deref(dist)
        ## cdef map[string,vector[int]].iterator it = m.begin()
        ## while it != m.end():
        ##     print deref(it).first.c_str()#, deref(it).second
        ##     inc(it)
        self.ptr.set_tip_conditionals(m)

    def optimize_global_dispersal_extinction(self, bool marginal, RateModel m):
        cdef double initL = super2double(self.ptr.eval_likelihood(marginal))
        print "Init L:", initL
        cdef _OptimizeBioGeo* opt = new _OptimizeBioGeo(self.ptr, m.ptr, marginal)
        cdef vector[double] disext = opt.optimize_global_dispersal_extinction()
        print "D:", disext[0]
        print "E:", disext[1]
        m.setup_D(disext[0])
        m.setup_E(disext[1])
        m.setup_Q()
        self.ptr.update_default_model(m.ptr)
        self.ptr.set_store_p_matrices(True)
        cdef double finalL = super2double(self.ptr.eval_likelihood(marginal))
        print "final L:", finalL
        self.ptr.set_store_p_matrices(False)
        
    def ancsplits(self, Tree intree, bool marginal, RateModel m, list areas):
        cdef int n = intree.ptr.getInternalNodeCount()
        cdef int i
        cdef _Node* node
        cdef map[vector[int],vector[_AncSplit]] ras
        cdef map[int,string] area_i2s
        cdef _BioGeoTreeTools tt

        self.ptr.set_use_stored_matrices(True)
        self.ptr.prepare_ancstate_reverse()

        for i, a in enumerate(areas):
            area_i2s[i] = string(<char *>a)

        for i in range(n):
            print i
            node = intree.ptr.getInternalNode(i)
            ras = self.ptr.calculate_ancsplit_reverse(deref(node), marginal)
            tt.summarizeSplits(node, ras, area_i2s, m.ptr)

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
    def __dealloc__(self):
        del self.ptr

    def read_treefile(self, filename):
        cdef string s = string(<char *>filename)
        cdef vector[_Tree*] *trees = new vector[_Tree*]() #deref(self.trees)
        cdef vector[_Tree*] tv
        self.ptr.readMultipleTreeFile(<string>s, deref(trees))
        tv = deref(trees)
        print tv.size()
