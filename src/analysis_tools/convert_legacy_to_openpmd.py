#!/usr/bin/env python

import argparse
import traceback
import synergia
import h5py


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a particle bunch stored in the legacy file format to a newer OpenPMD file"
    )
    parser.add_argument("-i", "--input", help="Input legacy bunch filename", type=str)
    parser.add_argument(
        "-o",
        "--output",
        help="Output OpenPMD bunch filename",
        type=str,
        default="bunch_openpmd.h5",
    )
    args = parser.parse_args()
    return args


def convert(filename_in: str, filename_out: str):
    file = h5py.File(filename_in)
    parts = file["particles"]
    size = parts.shape[0]

    # these will be reset upon read, but we need a
    # reference particle to create a bunch
    ref_part = synergia.foundation.Reference_particle(1, 1, 1)
    bunch = synergia.bunch.Bunch(ref_part, size, size)
    bunch.read_file_legacy(filename_in)

    bunch.write_openpmd_file(filename_out)
    return


if __name__ == "__main__":
    try:
        inputs = parse_args()
        convert(inputs.input, inputs.output)
    except Exception as e:
        print("Exception in user code:")
        traceback.print_exception(e)
