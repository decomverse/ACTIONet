from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools
import os
from pathlib import Path
import distutils.ccompiler
import sysconfig


__version__ = '1.0'

source_files = list()
for (dirpath, dirnames, filenames) in os.walk('src'):
	for file in filenames:
		if file.endswith(".cc") or file.endswith(".c"):
			source_files += [os.path.join(dirpath, file)]
source_files+=["python_interface/wrapper.cc"]



os.environ['CC'] = "ccache gcc"
#os.environ['LDFLAGS'] = "-Wl,--start-group " + os.environ['MKLROOT']+"/lib/intel64/libmkl_intel_ilp64.a " + os.environ['MKLROOT']+ "/lib/intel64/libmkl_sequential.a " + os.environ['MKLROOT'] + "/lib/intel64/libmkl_core.a" + " -Wl,--end-group"
#MKL_files=[os.environ['MKLROOT']+"/lib/intel64/libmkl_intel_ilp64.a", os.environ['MKLROOT']+ "/lib/intel64/libmkl_sequential.a", os.environ['MKLROOT'] + "/lib/intel64/libmkl_core.a "]

MKL=[os.environ['MKLROOT']+"/lib/intel64/libmkl_intel_ilp64.a", os.environ['MKLROOT']+ "/lib/intel64/libmkl_sequential.a", os.environ['MKLROOT'] + "/lib/intel64/libmkl_core.a"]


#extra_compile_args = sysconfig.get_config_var('CFLAGS').split()
extra_compile_args = ["-w", "-m64"]


def parallelCCompile(self, sources, output_dir=None, macros=None, include_dirs=None, debug=0, extra_preargs=None, extra_postargs=None, depends=None):
    # those lines are copied from distutils.ccompiler.CCompiler directly
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)
    # parallel code
    N=8 # number of parallel compilations
    import multiprocessing.pool
    def _single_compile(obj):
        try: src, ext = build[obj]
        except KeyError: return
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)
    # convert to list, imap is evaluated on-demand
    list(multiprocessing.pool.ThreadPool(N).imap(_single_compile,objects))
    return objects
distutils.ccompiler.CCompiler.compile=parallelCCompile


class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


ACTIONet_module = Extension(
        '_ACTIONet',
        source_files,
	define_macros=[('MKL_ILP64', None)],
        include_dirs=['include', os.environ['MKLROOT']+'/include', 'include/arma','include/ACTIONet', 'include/ACTIONet/hnsw', 'include/ACTIONet/leiden', 'include/ACTIONet/leiden/igraph','include/ACTIONet/HDBSCAN', 'include/ACTIONet/s_gd2', 'include/ACTIONet/uwot', get_pybind_include(), get_pybind_include(user=True)],
	library_dirs = [os.environ['MKLROOT']+'/lib/intel64'],
	libraries = ["pthread", "m", "dl"],
	extra_objects=MKL+MKL+MKL,
	#extra_compile_args = extra_compile_args,
	#extra_link_args = MKL_files,
        language='c++'
    )


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14/17] compiler flag.

    The newer version is prefered over c++11 (when it is available).
    """
    flags = ['-std=c++17', '-std=c++14', '-std=c++11']

    for flag in flags:
        if has_flag(compiler, flag): return flag

    raise RuntimeError('Unsupported compiler -- at least C++11 support '
                       'is needed!')


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
