from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# from Cython.Build import cythonize

ext = Extension(
    "lgcpp",
    ["lgcpp.pyx",],
    language = "c++",
    ## include_dirs = ["/home/rree/mypy/include"],
    ## library_dirs = ["/home/rree/mypy/lib"],
    libraries = ["lapack", "blas", "gfortran", "gsl", "gslcblas", "m",
                 "pthread", "nlopt", "armadillo"],
    extra_link_args = ['-fopenmp','-O2','-pthread','-fno-strict-aliasing','-g',
              '-Wall','-fPIC','-fwrapv'],
    extra_objects = [
        "AncSplit.o", 
        "InputReader.o", 
        "BayesianBioGeoAllDispersal.o", 
        "BayesianBioGeo.o", 
        "BioGeoTree.o", 
        "BioGeoTreeTools.o", 
        "BranchSegment.o", 
        "OptimizeBioGeo.o", 
        "OptimizeBioGeoAllDispersal.o", 
        "OptimizeBioGeoAllDispersal_nlopt.o", 
        "RateMatrixUtils.o", 
        "RateModel.o", 
        "Utils.o", 
        "node.o", 
        "tree.o", 
        "tree_reader.o", 
        "tree_utils.o", 
        "superdouble.o",
        "clock.o",
        "mataid.o",
        ## "blas.o",
        ## "lapack.o",
        "my_expokit.o",
        "my_matexp.o"
        ]
    )

## ext = Extension(
##     "lgcpp",
##     ["lgcpp.cpp",
##      "AncSplit.cpp", 
##      "InputReader.cpp", 
##      "BayesianBioGeoAllDispersal.cpp", 
##      "BayesianBioGeo.cpp", 
##      "BioGeoTree.cpp", 
##      "BioGeoTreeTools.cpp", 
##      "BranchSegment.cpp", 
##      "OptimizeBioGeo.cpp", 
##      "OptimizeBioGeoAllDispersal.cpp", 
##      "OptimizeBioGeoAllDispersal_nlopt.cpp", 
##      "RateMatrixUtils.cpp", 
##      "RateModel.cpp", 
##      "Utils.cpp", 
##      "node.cpp", 
##      "tree.cpp", 
##      "tree_reader.cpp", 
##      "tree_utils.cpp", 
##      "superdouble.cpp"],
##     language = "c++",
##     include_dirs = ["/home/rree/mypy/include"],
##     library_dirs = ["/home/rree/mypy/lib"],
##     libraries = "lapack blas gfortran gsl gslcblas m pthread nlopt python2.7".split(),
##     extra_link_args = ['-fopenmp','-O2','-pthread','-fno-strict-aliasing','-g',
##               '-Wall','-fPIC','-fwrapv'],
##     extra_objects = [
##         "clock.o",
##         "mataid.o",
##         "blas.o",
##         "lapack.o",
##         "my_expokit.o",
##         "my_matexp.o"
##         ]
##     )


if __name__ == "__main__":
    setup(name = "Lagrange C++ extension module",
          cmdclass = {"build_ext": build_ext},
          ext_modules = [ext])
      
