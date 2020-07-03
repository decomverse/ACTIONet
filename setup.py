from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools
import os
from pathlib import Path
import distutils.ccompiler
import sysconfig
from numpy.distutils.system_info import get_info
import pybind11
import multiprocessing

__version__ = '1.0'

N=multiprocessing.cpu_count()

def parallelCCompile(self, sources, output_dir=None, macros=None, include_dirs=None, debug=0, extra_preargs=None, extra_postargs=None, depends=None):
    # those lines are copied from distutils.ccompiler.CCompiler directly
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)
    # parallel code
    import multiprocessing.pool
    def _single_compile(obj):
        try: src, ext = build[obj]
        except KeyError: return
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)
    # convert to list, imap is evaluated on-demand
    list(multiprocessing.pool.ThreadPool(N).imap(_single_compile,objects))
    return objects
distutils.ccompiler.CCompiler.compile=parallelCCompile


def has_flag(compiler, flagname):
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True

def cpp_flag(compiler):
    flags = ['-std=c++17', '-std=c++14', '-std=c++11']

    for flag in flags:
        if has_flag(compiler, flag): return flag

    raise RuntimeError('Unsupported compiler -- at least C++11 support '
                       'is needed!')

                       
blas_opt = get_info('lapack_opt')

ACTIONet_source_files=["python_interface/wrapper.cc"]
for (dirpath, dirnames, filenames) in os.walk('src'):
    for file in filenames:
        if file.endswith(".cc") or file.endswith(".c"):
            ACTIONet_source_files+=[os.path.join(dirpath, file)]

ACTIONet_header_dirs=['include']
for (dirpath, dirnames, filenames) in os.walk('include'):
    for dirname in dirnames: ACTIONet_header_dirs += [os.path.join(dirpath, dirname)]

ACTIONet_macros=list()
if ("define_macros" in blas_opt): ACTIONet_macros+=blas_opt['define_macros']

ACTIONet_lib_dirs=list()
if ("library_dirs" in blas_opt): ACTIONet_lib_dirs+=blas_opt['library_dirs']

ACTIONet_libs=list()
if ("libraries" in blas_opt): ACTIONet_libs+=blas_opt['libraries']

ACTIONet_extra_args=list()
if ("extra_compile_args" in blas_opt): ACTIONet_extra_args+=blas_opt['extra_compile_args']

ACTIONet_extra_link_args=list()
if ("extra_link_args" in blas_opt): ACTIONet_extra_link_args+=blas_opt['extra_link_args']


os.environ['CC'] = "gcc"
os.environ['CFLAGS'] = " -w"

ACTIONet_module = Extension(
    '_ACTIONet',
    ACTIONet_source_files,
    define_macros=ACTIONet_macros,
    include_dirs=[pybind11.get_include(False), pybind11.get_include(True)]+ACTIONet_header_dirs,    
    library_dirs=ACTIONet_lib_dirs,
    libraries=ACTIONet_libs,
    extra_compile_args=ACTIONet_extra_args,
    extra_link_args=ACTIONet_extra_link_args,
    language='c++'
)


setup(
    name='ACTIONet',
    version=__version__,
    author='Shahin Mohammadi',
    author_email='shahin.mohammadi@gmail.com',
    url='https://github.com/shmohammadi86/ACTIONet',
    description='ACTIONet single-cell analysis framework',
    long_description='',
    packages=setuptools.find_packages("lib"),
    package_dir={"": "lib"},
    ext_modules=[ACTIONet_module],
    install_requires=['pybind11>=2.4', 'numpy'],
    setup_requires=['pybind11>=2.4', 'numpy'],
    zip_safe=False
)
