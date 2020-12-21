import multiprocessing
import os
import platform
import re
import sys
import subprocess
import distutils

from setuptools import find_packages, setup, Extension
from setuptools.command.build_ext import build_ext

__version__ = '1.0'




def parallelCCompile(self, sources, output_dir=None, macros=None, include_dirs=None, debug=0, extra_preargs=None, extra_postargs=None, depends=None):
    # those lines are copied from distutils.ccompiler.CCompiler directly
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)
    # parallel code
    N=multiprocessing.cpu_count() - 2
    import multiprocessing.pool
    def _single_compile(obj):
        try: src, ext = build[obj]
        except KeyError: return
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)
    # convert to list, imap is evaluated on-demand
    list(multiprocessing.pool.ThreadPool(N).imap(_single_compile,objects))
    return objects
distutils.ccompiler.CCompiler.compile=parallelCCompile

def read(path):
    with open(path, 'r') as f:
        return f.read()

def get_numpy_include():
    import numpy
    return numpy.get_include()


# Classes below were taken from https://github.com/pybind/cmake_example
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
            '-DPYTHON_EXECUTABLE=' + sys.executable,
            '-DNUMPY_INCLUDE_DIRS=' + get_numpy_include(),
        ]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j' + str(multiprocessing.cpu_count())]

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)


setup(
    name='ACTIONet',
    version=__version__,
    author='Shahin Mohammadi',
    author_email='shahin.mohammadi@gmail.com',
    url='https://github.com/shmohammadi86/ACTIONet',
    description='ACTIONet single-cell analysis framework',
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    keywords='',
    python_requires='>=3.5',
    packages=find_packages(),
    ext_modules=[CMakeExtension('_ACTIONet')],
    cmdclass=dict(build_ext=CMakeBuild),
    install_requires=read('requirements.txt').strip().split('\n'),
    setup_requires=['numpy'],
    zip_safe=False,
    include_package_data=True,
    classifiers=[],
)
