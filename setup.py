import distutils.sysconfig as sysconfig
import multiprocessing
import os
import platform
import re
import subprocess
import sys
from distutils.sysconfig import get_python_inc
from distutils.version import LooseVersion

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext

from setup_helpers import ParallelCompile, naive_recompile

ParallelCompile("NPY_NUM_BUILD_JOBS", needs_recompile=naive_recompile).install()


__version__ = "0.4.1"


def read(path):
    with open(path, "r") as f:
        return f.read()


def get_numpy_include():
    import numpy

    return numpy.get_include()


# Classes below were taken from https://github.com/pybind/cmake_example
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " + ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r"version\s*([\d.]+)", out.decode()).group(1))
            if cmake_version < "3.1.0":
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep
        cmake_args = ["-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir, "-DPYTHON_EXECUTABLE=" + sys.executable, "-DNUMPY_INCLUDE_DIRS=" + get_numpy_include(), "-DPYTHON_INCLUDE_DIR=" + get_python_inc(), "-DPYTHON_LIBRARY=" + sysconfig.get_config_var("LIBDIR")]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            cmake_args += ["-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            build_args += ["--", "-j" + str(multiprocessing.cpu_count())]

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get("CXXFLAGS", ""), self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=self.build_temp)


setup(
    name="ACTIONet",
    version=__version__,
    author="Shahin Mohammadi",
    author_email="shahin.mohammadi@gmail.com",
    url="https://github.com/shmohammadi86/ACTIONet",
    description="ACTIONet single-cell analysis framework",
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    keywords="",
    python_requires=">=3.8",
    packages=find_packages(),
    ext_modules=[CMakeExtension("_ACTIONet")],
    cmdclass=dict(build_ext=CMakeBuild),
    install_requires=read("requirements.txt").strip().split("\n"),
    setup_requires=["numpy"],
    zip_safe=False,
    include_package_data=True,
    classifiers=[],
    options={"bdist_wheel":{"universal": True}}
)
