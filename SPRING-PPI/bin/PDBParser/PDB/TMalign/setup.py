from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize(
        "tmalign.pyx",
        sources = ["TMalign_wrapper.cpp"],
        language = "c++",
        extra_compile_args = ["-static", "-O3" , "-ffast-math","-lm"],
    ))
