#!/usr/bin/env python

from macro_bunch import *
from macro_bunch_store import foo
import unittest
class Test_Macro_bunch(unittest.TestCase):
    def test_01_construct(self):
        mb = Macro_bunch()

    def test_02_init_test(self):
        mb = Macro_bunch()
        mb.init_test()

    def test_03(self):
        foo(1)
        foo(0)

if __name__ == '__main__':
    macro_bunch_suite = unittest.TestLoader().loadTestsFromTestCase(Test_Macro_bunch)
    unittest.TextTestRunner(verbosity=2).run(macro_bunch_suite)
