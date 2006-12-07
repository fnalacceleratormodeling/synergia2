#!/usr/bin/env python

from macro_bunch import *
import unittest
class Test_Macro_bunch(unittest.TestCase):
    def test_01_construct1(self):
        mb = Macro_bunch()

if __name__ == '__main__':
    macro_bunch_suite = unittest.TestLoader().loadTestsFromTestCase(Test_Macro_bunch)
    unittest.TextTestRunner(verbosity=2).run(macro_bunch_suite)
