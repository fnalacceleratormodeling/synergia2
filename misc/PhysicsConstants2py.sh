#!/bin/sh

fnal_dir=/home3/amundson/work/fnal/branches/jfa1
grep '#define.*[0-9]' $fnal_dir/basic_toolkit/include/PhysicsConstants.h \
| sed 's/#define \+\([^ ]*\) \+/\1=/' \
| sed 's%/\*\(.*\)\*/%#\1%' \
> physics_constants.py
