from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension("graph_python",
              sources=['graph_python.pyx'],
              language='c++',
              extra_compile_args=['-std=c++17'])
]

setup(name='Test', ext_modules=cythonize(extensions, compiler_directives={'language_level': "3"}))
