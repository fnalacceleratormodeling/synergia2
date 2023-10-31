#!/usr/bin/env python

import argparse
import traceback
import numpy as np
from dataclasses import dataclass
from enum import Enum, unique
import openpmd_api as io

# constant, speed of light
C = 299792458.0


@unique
class Coords(Enum):
    x = str("position_x")
    xp = str("moments_x")
    y = str("position_y")
    yp = str("moments_y")
    cdt = str("position_z")
    dpop = str("moments_z")
    pz = str("pz")
    energy = str("energy")
    t = str("t")
    z = str("z")

    # thanks to https://stackoverflow.com/a/55500795
    # magic methods for argparse compatibility
    def __str__(self):
        return self.name.lower()

    def __repr__(self):
        return str(self)

    @staticmethod
    def argparse(s):
        try:
            return Coords[s]
        except KeyError:
            return s


@dataclass
class Options:
    hcoord: Coords
    vcoord: Coords
    inputfile: str
    iteration: int
    num_bins: int
    minh: float
    maxh: float
    minv: float
    maxv: float
    outputfile: str = ""


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="OpenPMD-series-filename", type=str)

    parser.add_argument(
        "--iteration", help="iteration of OpenPMD series to plot", type=int, default=0
    )
    parser.add_argument(
        "--minh", help="minimum limit on horizontal axis data", type=float, 
        default=-1.0 * np.finfo(np.float64).max
    )
    parser.add_argument(
        "--maxh", help="maximum limit on horizontal axis data", type=float,
        default=np.finfo(np.float64).max
    )
    parser.add_argument(
        "--minv", help="minimum limit on vertical axis data", type=float,
        default=-1.0 * np.finfo(np.float64).max
    )
    parser.add_argument(
        "--maxv", help="maximum limit on vertical axis data", type=float,
        default=np.finfo(np.float64).max

    )
    parser.add_argument(
        "--bins",
        help="number of bins for histogram to draw",
        default=25,
        type=int,
    )
    parser.add_argument(
        "--output", help="save output to file (not on by default)", type=str
    )
    parser.add_argument(
        "xcoord", help="x-coord to plot", type=Coords.argparse, choices=list(Coords)
    )
    parser.add_argument(
        "ycoord", help="y-coord to plot", type=Coords.argparse, choices=list(Coords)
    )

    args = parser.parse_args()
    return args


def do_plots(opts: Options):
    import seaborn as sns
    import matplotlib.pyplot as plt

    sns.set_theme(style="darkgrid")

    series = io.Series(opts.inputfile, io.Access_Type.read_only)
    i = series.iterations[opts.iteration]
    parts = i.particles["bunch_particles"]
    parts_df = parts.to_df()

    mass = parts.get_attribute("mass")
    beta = parts.get_attribute("beta_ref")

    print(
        f'''-------- Using the following parameters ----------
iteration: {opts.iteration}; xcoord: {opts.hcoord.value}; ycoord: {opts.vcoord.value};
minh, maxh: {opts.minh}/{opts.maxh};
minv, maxv: {opts.minv}/{opts.maxv};
num-bins: {opts.num_bins}
--------------------------------------------------'''
    )

    if opts.hcoord == Coords("pz") or opts.vcoord == Coords("pz"):
        parts_df["pz"] = p_ref + parts_df["moments_z"]

    if opts.hcoord == Coords("energy") or opts.vcoord == Coords("energy"):
        pz = p_ref + parts_df["moments_z"]
        parts_df["energy"] = np.sqrt(pz * pz + mass**2)

    if opts.hcoord == Coords("t") or opts.vcoord == Coords("t"):
        parts_df["t"] = parts_df["position_z"] * 1.0e9 / C

    if opts.hcoord == Coords("z") or opts.vcoord == Coords("z"):
        parts_df["z"] = parts_df["position_z"] * beta

    to_plot_df = parts_df[
        (parts_df[opts.hcoord.value] > opts.minh)
        & (parts_df[opts.hcoord.value] < opts.maxh)
        & (parts_df[opts.vcoord.value] > opts.minv)
        & (parts_df[opts.vcoord.value] < opts.maxv)
    ]

    if to_plot_df.empty:
        raise Exception("Empty data frame to plot after applying co-ordinate selection."+"\n"+"Please check input co-ordinates and run again!")

    g = sns.JointGrid(
        data=to_plot_df, x=opts.hcoord.value, y=opts.vcoord.value, marginal_ticks=True
    )
    g.plot_joint(sns.histplot, thresh=None, cmap="mako", bins=opts.num_bins)
    g.plot_marginals(sns.histplot, kde=False)

    plt.suptitle("Synergia3 Phase Space Distribution", fontsize="medium", y=0.985)
    if opts.outputfile:
        g.savefig(opts.outputfile)
    else:
        plt.show()

    return


if __name__ == "__main__":
    try:
        inputs = parse_args()
        opts = Options(
            hcoord=inputs.xcoord,
            vcoord=inputs.ycoord,
            minh=inputs.minh,
            maxh=inputs.maxh,
            minv=inputs.minv,
            maxv=inputs.maxv,
            inputfile=inputs.filename,
            iteration=inputs.iteration,
            num_bins=inputs.bins,
            outputfile=inputs.output,
        )
        do_plots(opts)
    except Exception as e:
        print("Exception in user code:")
        traceback.print_exception(e)
