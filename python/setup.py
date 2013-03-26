# Author: Vlad Morariu
# Created: 2012-01-25
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
import sys

if not sys.platform == 'win32':
   include_dirs=['../include', numpy.get_include()]
   library_dirs=['../lib']
   libraries=['figtree', 'ann_figtree_version']
   extra_compile_args=['-DNDEBUG', '-O2']
else:
   include_dirs=['../include', numpy.get_include()]
   library_dirs=['../lib']
   libraries=['figtree', 'ann_figtree_version']
   extra_compile_args=['-O2']

ext_modules = [Extension(
                         "figtree",
                         ["figtree.pyx"],
                         language="c++",
                         include_dirs=include_dirs,
                         library_dirs=library_dirs,
                         libraries=libraries,
                         extra_compile_args=extra_compile_args,
                         extra_link_args=[])]

for e in ext_modules:
   e.pyrex_directives={'embedsignature': True}
   
setup(
    ext_modules=ext_modules,
    cmdclass={'build_ext' : build_ext}
)
