CHEF_ROOT="/home/egstern/contract-chef"
INSTALLDIR="${CHEF_ROOT}/install"
BOOST_INC=${INSTALLDIR}/include/boost-1_34_1
CHEF_BIN=${INSTALLDIR}/bin
CHEF_LIB=${INSTALLDIR}/lib
if ! echo ${PATH} | grep -q ${CHEF_BIN} ; then
    PATH=${CHEF_BIN}:${PATH}
fi

export CHEF_ROOT
export INSTALLDIR
export BOOST_INC
export CHEF_BIN
export CHEF_LIB

if [ -z ${LD_LIBRARY_PATH} ] ; then
    export LD_LIBRARY_PATH=${CHEF_LIB}
else
    if ! echo ${LD_LIBRARY_PATH} | grep -q ${CHEF_LIB} ; then
	export LD_LIBRARY_PATH=$CHEF_LIB}:${LD_LIBRARY_PATH}
    fi
fi
