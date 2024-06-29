#!/usr/bin/env python
import pytest
import synergia
from synergia.lattice import Lattice, Lattice_element
from synergia.foundation import Reference_particle

name = "foo"

def test_evil_functon():
    lattice = Lattice(name)
    try:
        lattice.do_not_call_this_function()
    except:
        raise Runtime_Error('I told you not to call this function!!!')
    
    assert True
