#!/usr/bin/env python

import argparse
from dataclasses import dataclass
from enum import Enum, unique
from typing import List
import traceback
import openpmd_api as io
import numpy as np
from cycler import Cycler, cycler
import matplotlib.pyplot as plt


@unique
class Coords(Enum):
    x = [str("position"), str("x")]
    xp = [str("moments"), str("x")]
    y = [str("position"), str("y")]
    yp = [str("moments"), str("y")]
    z = [str("position"), str("z")]
    zp = [str("moments"), str("z")]
    pz = [str("moments"), str("z"), str("pz")]
    energy = [str("moments"), str("z"), str("energy")]

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
    indices: List[int]
    coords: List[Coords]
    inputfile: str
    userep: bool = False
    oneplot: bool = False
    outputfile: str = ""


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="OpenPMD-series-filename", type=str)

    parser.add_argument(
        "--indices", nargs="*", help="indices to plot", type=int, default=[0]
    )

    parser.add_argument(
        "--output", help="save output to file (not on by default)", type=str
    )
    parser.add_argument(
        "--oneplot", help="put all plots on the same axis", action="store_true"
    )

    parser.add_argument(
        "--coords",
        nargs="+",
        help="coord to plot",
        type=Coords.argparse,
        choices=list(Coords),
        action="extend",
    )

    args = parser.parse_args()
    return args


def get_layout(num):
    if num == 1:
        return 1, 1
    elif num == 2:
        return 2, 1
    elif num == 3:
        return 3, 1
    elif num == 4:
        return 2, 2
    elif num <= 6:
        return 3, 2
    elif num <= 9:
        return 3, 3
    elif num <= 12:
        return 3, 4
    elif num <= 16:
        return 4, 4
    else:
        raise ValueError("Too many plots!")


def do_plot_coord(
    filename: str,
    coord: Coords,
    indices: List[int],
    line_style: Cycler,
    ax: plt.Axes,
):
    series = io.Series(filename, io.Access_Type.read_only)
    p_ref: float = series.get_attribute("pz")
    mass: float = series.get_attribute("mass")

    dim = len(series.iterations)
    x = np.zeros(dim)
    y = np.zeros(dim)
    _xlabel: str = ""
    for count, iteration in series.iterations.items():
        x[count] = iteration.get_attribute("track_s")
        _xlabel = "track_s"
    del series

    color_cycle = cycler(color=["c", "m", "y", "k"])

    for color_style, index in zip(color_cycle, indices):
        labelstr: str = ""
        series = io.Series(filename, io.Access_Type.read_only)

        for count, iteration in series.iterations.items():
            _parts = iteration.particles["track_coords"]
            _data = _parts[coord.value[0]][coord.value[1]]
            y[count] = _data[index]
            series.flush()
            if coord == Coords.pz:
                y = (y + 1) * p_ref
            if coord == Coords.energy:
                y = (y + 1) * p_ref
                y = np.sqrt(y * y + mass * mass)
        labelstr = str(coord) + str(" ") + str(index)
        _style = {**color_style, **line_style}
        ax.plot(x, y, **_style, label=labelstr)
        ax.set_xticks(x)
        ax.set_xlabel(_xlabel)

    return


def do_plots(opts: Options):
    num_plots = len(opts.coords)
    rows, cols = get_layout(num_plots)
    plt.rc("lines", linewidth=2)
    style_cycle: Cycler = (
        cycler(linestyle=["solid", "dashed", "dashdot", "dotted"])
    ) + (cycler(marker=["o", "*", "+", "^"]))

    if opts.oneplot:
        fig, ax = plt.subplots()
        fig.suptitle("Synergia3 Track Viewer", fontsize="medium")
    else:
        fig, ax = plt.subplots(rows, cols)
        fig.suptitle("Synergia3 Track Viewer", fontsize="medium")
        if (not isinstance(ax, plt.Axes)) and len(ax) != 1:
            ax = ax.flatten()
    if isinstance(ax, plt.Axes):
        for coord_to_plot, linestyle in zip(opts.coords, style_cycle):
            do_plot_coord(opts.inputfile, coord_to_plot, opts.indices, linestyle, ax)
    else:
        for coord_to_plot, linestyle, _ax in zip(opts.coords, style_cycle, ax):
            do_plot_coord(opts.inputfile, coord_to_plot, opts.indices, linestyle, _ax)

    if isinstance(ax, plt.Axes):
        ax.legend()
    else:
        for _ax in ax:
            _ax.legend()

    if opts.outputfile:
        plt.savefig(opts.outputfile)
    else:
        plt.show()

    return


if __name__ == "__main__":
    try:
        inputs = parse_args()
        opts = Options(
            indices=inputs.indices,
            inputfile=inputs.filename,
            oneplot=inputs.oneplot,
            outputfile=inputs.output,
            coords=inputs.coords,
        )
        do_plots(opts)

    except Exception as e:
        print("Exception in user code:")
        traceback.print_exception(e)
