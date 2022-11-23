#!/usr/bin/env python
import sys
sys.path.append('../../..')
import local_paths

from mpi4py import MPI
import synergia
from synergia.utils import Logger
from nose.tools import *

def test_construct1():
    logger = Logger(0)

def test_construct2():
    logger = Logger(0,"logger-log-py-2")

def test_construct3():
    logger = Logger("logger-log-py-3")

def test_doit1():
    Logger(0).write("doit1\n")

def test_doit2():
    Logger(0,"logger-log-py-2").write("doit2\n")

def test_doit1():
    Logger("logger-log-py-3").write("doit3\n")
