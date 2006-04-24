#!/usr/bin/env python
import os.path
gourmet_dir = "/home2/amundson/work/contract-synergia2/build/chef/python-bindings/src/.libs"
impact_dir = "/home2/amundson/work/contract-synergia2/build/impact/Forthon_Interfaces"
synergia2_dir = "/home2/amundson/work/contract-synergia2/build/synergia2"
octapy_dir = os.path.join(synergia2_dir,"octapy")

import sys
sys.path.append(gourmet_dir)
sys.path.append(impact_dir)
sys.path.append(synergia2_dir)
sys.path.append(octapy_dir)

