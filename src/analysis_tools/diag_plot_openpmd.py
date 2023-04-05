#!/usr/bin/env python

import argparse
from dataclasses import dataclass
from cycler import Cycler, cycler
from enum import Enum, unique
from itertools import product, combinations
from typing import List, Tuple, Union
import traceback
import openpmd_api as io
import numpy as np
import matplotlib.pyplot as plt

idx_dict = {"x": 0, "xp": 1, "y": 2, "yp": 3, "z": 4, "zp": 5}
diags: List[Tuple[str, Union[str, Tuple[str, int]]]] = [
    (str("x_emit"), str("emitx")),
    (str("y_emit"), str("emity")),
    (str("z_emit"), str("emitz")),
    (str("xy_emit"), str("emitxy")),
    (str("xyz_emit"), str("emitxyz")),
    (str("particles"), str("num_particles")),
]
for p in product(["x", "xp", "y", "yp", "z", "zp"], ["std", "mean"]):
    state = "std" if p[1] == "std" else "mean"
    diags.append((str(p[0] + "_" + p[1]), (str(state), idx_dict[p[0]])))

for c in combinations(["x", "xp", "y", "yp", "z", "zp"], 2):
    for p in product([str(c[0] + "_" + c[1] + "_")], ["corr", "mom2"]):
        state = "corr" if p[1] == "corr" else "mom2"
        diags.append(
            (str(p[0] + p[1]), (str(state), idx_dict[c[0]] * 6 + idx_dict[c[1]]))
        )


# thanks to https://stackoverflow.com/a/70294405
def inject_items(d, items):
    for k, v in items:
        d[k] = v


@unique
class Diags(Enum):
    inject_items(locals(), diags)

    # thanks to https://stackoverflow.com/a/55500795
    # magic methods for argparse compatibility
    def __str__(self):
        return self.name.lower()

    def __repr__(self):
        return str(self.name)

    @staticmethod
    def argparse(s):
        try:
            return Diags[s]
        except KeyError:
            return s


@dataclass
class Options:
    diags_to_plot: List[Diags]
    inputfiles: List[str]
    userep: bool = False
    oneplot: bool = False
    outputfile: str = ""


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="OpenPMD-series-filename", type=str)

    parser.add_argument(
        "--output", help="save output to file (not on by default)", type=str
    )
    parser.add_argument(
        "--userep",
        help="use repetition instead of s for independent variable",
        action="store_true",
    )
    parser.add_argument(
        "--oneplot", help="put all plots on the same axis", action="store_true"
    )

    parser.add_argument(
        "diag",
        nargs="*",
        help="diagnostic to plot",
        type=Diags.argparse,
        choices=list(Diags),
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


def do_plot_diag(
    diag: Diags,
    filenames: List[str],
    line_style: Cycler,
    _ax: plt.Axes,
):
    color_cycle = cycler(color=["c", "m", "y", "k"])
    for color_style, filename in zip(color_cycle, filenames):
        series = io.Series(filename, io.Access_Type.read_only)
        dim = len(series.iterations)
        x = np.zeros(dim)
        y = np.zeros(dim)
        _xlabel: str = ""
        for count, iteration in series.iterations.items():
            if opts.userep:
                x[count] = iteration.get_attribute("repetition")
                _xlabel = "repetition"
            else:
                x[count] = iteration.get_attribute("s")
                _xlabel = "s"

        labelstr: str = ""
        if len(filenames) != 1:
            labelstr = labelstr + filename

        if isinstance(diag.value, str):
            labelstr = labelstr + " " + diag.value
        elif isinstance(diag.value, tuple) and len(diag.value) == 2:
            labelstr = labelstr + " " + diag.value[0]

        for count, iteration in series.iterations.items():
            if isinstance(diag.value, str):
                y[count] = iteration.get_attribute(diag.value)
            elif isinstance(diag.value, tuple) and len(diag.value) == 2:
                y[count] = iteration.get_attribute(diag.value[0])[diag.value[1]]

        _style = {**color_style, **line_style}
        _ax.plot(x, y, **_style, label=labelstr)
        _ax.set_xticks(x)
        _ax.set_xlabel(_xlabel)
        
    return

def do_plots(opts: Options):
    num_plots = len(opts.diags_to_plot)
    rows, cols = get_layout(num_plots)
    plt.rc("lines", linewidth=2)
    style_cycle: Cycler = (
        cycler(linestyle=["solid", "dashed", "dashdot", "dotted"])
    ) + (cycler(marker=["o", "*", "+", "^"]))

    if opts.oneplot:
        fig, ax = plt.subplots()
        fig.suptitle("Synergia3 Phase Space Distribution", fontsize="medium")
    else:
        fig = plt.figure()
        fig.suptitle("Synergia3 Phase Space Distribution", fontsize="medium")
        gs = fig.add_gridspec(nrows=rows, ncols=cols, hspace=0)
        ax = gs.subplots(sharex="col")
        if (not isinstance(ax, plt.Axes)) and len(ax) != 1:
            ax = ax.flatten()
    if isinstance(ax, plt.Axes):
        for diag, linestyle in zip(opts.diags_to_plot, style_cycle):
            do_plot_diag(diag, opts.inputfiles, linestyle, ax)
    else:
        for diag, linestyle, _ax in zip(opts.diags_to_plot, style_cycle, ax):
            do_plot_diag(diag, opts.inputfiles, linestyle, _ax)

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
            diags_to_plot=inputs.diag,
            inputfiles=inputs.filename.split(","),
            userep=inputs.userep,
            oneplot=inputs.oneplot,
            outputfile=inputs.output,
        )
        do_plots(opts)

    except Exception as e:
        print("Exception in user code:")
        traceback.print_exception(e)
