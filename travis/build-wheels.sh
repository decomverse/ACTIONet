#!/bin/bash
yum -y install epel-release
yum config-manager --enable epel
yum -y install hdf5-devel
yum -y install suitesparse-devel nss blas-devel lapack-devel xorg-x11-server-Xvfb
python_versions=( cp38-cp38 cp39-cp39 cp310-cp310 )
cpython_versions=( cpython-3.8.14 cpython-3.9.14 cpython-3.10.7 )
python_major_versions=( 3.8 3.9 3.10 )
total=${#python_versions[*]}
for ((i=0; i<=$(( $total - 1 )); i++ ))
do
    cur_python=${python_versions[$i]}
    cur_cpython=${cpython_versions[$i]}
    cur_major_python=${python_major_versions[$i]}

    export PATH=/io/cmake:/opt/python/${cur_python}/bin:/usr/share/Modules/bin:/opt/rh/gcc-toolset-11/root/usr/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
    set -e -u -x

    function repair_wheel {
	wheel="$1"
	if ! auditwheel show "$wheel"; then
            echo "Skipping non-platform wheel $wheel"
	else
            auditwheel repair "$wheel" --plat "$PLAT" -w /io/wheelhouse/
	fi
    }


    # Install a system package required by our library
    yum install -y atlas-devel

    export PYTHON_INCLUDE_DIR=/opt/_internal/${cur_cpython}/include/python${cur_major_python}
    export PYTHON_LIBRARY=/opt/_internal/${cur_cpython}/lib
    # Compile wheels
    for PYBIN in /opt/python/${cur_python}/bin; do
	"${PYBIN}/pip" install -r /io/requirements.txt
	"${PYBIN}/pip" wheel /io/  -w wheelhouse/
    done

# Bundle external shared libraries into the wheels
    export LD_LIBRARY_PATH=/opt/_internal/${cur_cpython}/lib/python${cur_major_python}/site-packages:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=/opt/_internal/${cur_cpython}/lib/python${cur_major_python}/site-packages/numpy.libs:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=/opt/_internal/${cur_cpython}/lib/python${cur_major_python}/site-packages/scipy.libs:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=/opt/_internal/${cur_cpython}/lib/python${cur_major_python}/site-packages/Pillow.libs:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=/opt/_internal/${cur_cpython}/lib/python${cur_major_python}/site-packages/pyzmq.libs:$LD_LIBRARY_PATH
    for whl in wheelhouse/ACTIONet-0.3.0-${cur_python}-linux_x86_64.whl; do
	repair_wheel "$whl"
    done
    # Install packages and test
    for PYBIN in /opt/python/${cur_python}/bin/; do
	"${PYBIN}/pip" install ACTIONet --no-index -f /io/wheelhouse
    done
done
