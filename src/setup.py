from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext = Extension(
    "lgcpp",
    ["lgcpp.pyx",
     "AncSplit.cpp", 
     "InputReader.cpp", 
     "BayesianBioGeoAllDispersal.cpp", 
     "BayesianBioGeo.cpp", 
     "BioGeoTree.cpp", 
     "BioGeoTreeTools.cpp", 
     "BranchSegment.cpp", 
     "OptimizeBioGeo.cpp", 
     "OptimizeBioGeoAllDispersal.cpp", 
     "OptimizeBioGeoAllDispersal_nlopt.cpp", 
     "RateMatrixUtils.cpp", 
     "RateModel.cpp", 
     "Utils.cpp", 
     "node.cpp", 
     "tree.cpp", 
     "tree_reader.cpp", 
     "tree_utils.cpp", 
     "superdouble.cpp"],
    language = "c++",
    include_dirs = ["/home/rree/mypy/include"],
    library_dirs = ["/home/rree/mypy/lib"],
    libraries = "lapack blas gfortran gsl gslcblas m pthread nlopt".split(),
    extra_link_args = ["-fopenmp"],
    extra_objects = ["clock.o",
                     "my_expokit.o",
                     "mataid.o",
                     "blas.o",
                     "lapack.o",
                     "my_matexp.o"]
    )


if __name__ == "__main__":
    setup(name = "Lagrange C++ extension module",
          cmdclass = {"build_ext": build_ext},
          ext_modules = [ext])
      
