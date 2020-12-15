from distutils.core import setup, Extension
import numpy
setup(name = 'CustomConfLim', version = '1.0',  \
   ext_modules = [Extension('CustomConfLim', ['CustomConfLim.c'])],
   include_dirs=[numpy.get_include()])