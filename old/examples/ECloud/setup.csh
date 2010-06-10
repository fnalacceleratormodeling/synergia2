if (${?PATH}) then
    setenv PATH "/local/lebrun/Synergia/cca/install/bin:${PATH}"
else
    setenv PATH "/local/lebrun/Synergia/cca/install/bin"
endif

if (${?LD_LIBRARY_PATH}) then
    setenv LD_LIBRARY_PATH "/local/lebrun/Synergia/cca/install/lib:${LD_LIBRARY_PATH}"
else
    setenv LD_LIBRARY_PATH "/local/lebrun/Synergia/cca/install/lib"
endif
#
setenv PYTHONPATH ./
setenv PYTHONPATH ./ECloudCC:$PYTHONPATH
setenv PYTHONPATH /local/lebrun/Tech-Xlib/install/lib/python2.4/site-packages:$PYTHONPATH
setenv PYTHONPATH /local/lebrun/Tech-Xlib/install/lib:$PYTHONPATH
setenv PYTHONPATH /local/lebrun/Synergia/cca/install/lib:$PYTHONPATH
setenv PYTHONPATH /local/lebrun/Synergia/cca/install/lib/python2.4:$PYTHONPATH
setenv PYTHONPATH /local/lebrun/Synergia/cca/build/synergia2:$PYTHONPATH
setenv PYTHONPATH /local/lebrun/Synergia/cca/build/synergia2/s2_fish:$PYTHONPATH
setenv PYTHONPATH /local/lebrun/Synergia/cca/install/include/python2.4:$PYTHONPATH
#
setenv LD_LIBRARY_PATH /local/lebrun/Tech-Xlib/install/lib:${LD_LIBRARY_PATH}
#
setenv PATH /local/lebrun/Synergia/cca/install/bin:${PATH}
#
setenv HOMEELEC /local/lebrun/
#
