synergia module
==================

.. automodule:: synergia

.. autoclass:: Gourmet
    :members:
    :undoc-members:

.. autoclass:: Beam_parameters
    :members:
    :undoc-members:

.. autofunction:: propagate

.. autoclass:: Options
    :members:
    :undoc-members:

.. autoclass:: Diagnostics
    :members:
    :undoc-members:
    
    Diagnostics include (len index runs over number of steps, ix is the coordinate index:x=0,x'=1,y=2,y'=3,t=4,t'=5):
        s[len]
            total distance traveled by reference particle
        means[len,ix]
            mean beam position
        mom2s[len,ix1,ix2]
            second moment of beam distribution
        corrs[len,ix1,ix2]
            correlation coefficients of the beam distribution
        diagmom4s[len,ix]
            diagonal of the fourth moments of the beam distribution
        stds[len,ix]
            standard deviations of the beam distribution
        emitxs[len]
            horizontal emittance
        emitys[len]
            vertical emittance
        emitzs[len]
            longitudinal emittance
        emitxys[len]
            transverse (4d) emittance
        emitxyzs[len]
            full (6d) emittance
        u[6]
            unit conversion
        n[len]
            number of macro particles in bunch
        

.. autoclass:: Error_eater
    :members:
    :undoc-members:

    Eats output sent to stderr by C++ classes. Useful in conjunction with CHEF, which sends a great deal of extraneous output to stderr.

.. automodule:: synergia.matching
.. autofunction:: match_twiss_width
.. autofunction:: match_twiss_emittance
.. autofunction:: get_alpha_beta
.. autofunction:: envelope_match
