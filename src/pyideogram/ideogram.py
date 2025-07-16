#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 09:09:33 2023

@author: balthasar
"""
import matplotlib.path as mpath
import numpy as np
import lzma
import pickle
from operator import sub
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Patch
from matplotlib import cbook
from .matplotlib_extension import (
    Chromosomeconnection,
    SideRound,
    AngledBox,
    ReverseFancyPatch,
)


import warnings
from . import annotation
from .dataloader import load_cytobands
import matplotlib.patheffects as pe

import importlib.resources as pkg_resources


###############################################################################
# %% ################################ functions ################################
###############################################################################

BANDCOL = {
    "gneg": (0.95, 0.95, 0.95),
    "gpos25": (0.7, 0.7, 0.7),
    "gpos50": (0.5, 0.5, 0.5),
    "gpos75": (0.3, 0.3, 0.3),
    "gpos100": (0.1, 0.1, 0.1),
    "acen": (1, 0.38, 0.36),
    "gvar": (0.6, 0.6, 0.8),
    "stalk": (0.45, 0.4, 0.8),
}


def ideogramh(
    chrom: str,
    bands=None,
    ax: plt.Axes = None,
    round_pct: float = 0.55,
    color: dict = BANDCOL,
    style="round",
    names=False,
    textkwargs={},
    **ideokwargs
):
    """
    Plot the ideogram of a chromosome

    Parameters
    ----------
    chrom : str
        Mandatory. The chromosome which should be plotted.
    bands : None, path like, dictionary
        If this is None the cytoband file shipped with this package is used,
        which contains annotation for hg19 and hg39. They are differentiated by
        the absence or presence of 'chr' respectively.
        If bands is a dictionary, it is assumed to be loaded by the load_cytobands
        function from this package.
        If bands is a string or path and can be interpreted by open it load_cytobands
        will be executed to load the file inside the function instead.
    ax : plt.Axes
        Axes to plot on. If not given it is automatically determined.
    round_pct : float
        The rounding strength in pct of the rounding of the short side, which is 90%.
        This is scaled by the ratio and should look good at ~ 0.5
    color : dict
        A dictionary that maps the last column of a cytoband file to a color to
        display. Normally the last color contains some staining info, but you
        can put whatever you want and visualize anything combined with the color
        dict. You can get the default with pyideogram.BANDCOL
    style : str
        The style of the centromers. Currently 'round' and 'angled' are supported.
    names : bool
        Should the band name be printed in the center of the band. This is most
        likely only useful on moderate zoom levels into the ideogram. This is a
        bit WIP.
    textkwargs : dict
        A dictionary that can contain key word arguments for the matplotlib text
        method if names is True.
    barkwargs : dict
        A dictionary that can contain key word arguments for the matplotlib bar
        method.
    """

    # centromer style can be round or angled
    if bands is None:
        with pkg_resources.as_file(
            pkg_resources.files(annotation) / "cytoBands.pickle.xz"
        ) as f:
            bands = pickle.load(lzma.open(f))
    elif isinstance(bands, dict):
        pass
    else:
        try:
            bands = load_cytobands(bands)
        except:
            ValueError(
                """bands needs to be one of:
                       "None", "file like" or a dictionary that was loaded by the
                       "load_cytobands" function from this package"""
            )

    bands = bands[chrom]
    band_end = np.array(bands[0])
    band_start = np.array([0] + bands[0][:-1])
    band_name = bands[1]
    band_type = np.array(bands[2])
    band_num = len(band_end)

    if ax == None:
        ax = plt.gca()

    ideokwargs = cbook.normalize_kwargs(ideokwargs, Patch)
    ideokwargs = {"edgecolor": "k"} | ideokwargs
    # ideokwargs.setdefault('orientation', 'horizontal')
    orientation = ideokwargs.pop('orientation', "horizontal")

    # TODO: it should be possible to just hand orientation to ax.bar, as per the
    # undocumented api in the matplotlib source code but for some reason that
    # does not seem to work

    if orientation == "vertical":
        if ax.get_autoscale_on():
            ax.set_ylim(-0.01 * max(band_end), 1.01 * max(band_end))

        barpatches = ax.bar(
            [chrom] * band_num,
            band_end - band_start,
            1,
            band_start,
            # maybe defaultdict is better.
            color=[color[bt] if bt in color else (
                1, 1, 1) for bt in band_type],
            ** ideokwargs
        )

        cornerref = ([0, 1], [2, 3])

    else:
        if ax.get_autoscale_on():
            ax.set_xlim(-0.01 * max(band_end), 1.01 * max(band_end))

        barpatches = ax.barh(
            [chrom] * band_num,
            band_end - band_start,
            1,
            band_start,
            # maybe defaultdict is better.
            color=[color[bt] if bt in color else (
                1, 1, 1) for bt in band_type],
            ** ideokwargs
        )

        cornerref = ([0, 3], [1, 2])

    # remove all arguments that are used in bar internally.
    ideokwargs.pop('edgecolor', None)
    ideokwargs.pop('linewidth', None)
    ideokwargs.pop('hatch', None)
    ideokwargs.pop('xerr', None)
    ideokwargs.pop('yerr', None)
    ideokwargs.pop('error_kw', None)
    ideokwargs.pop('ecolor', None)
    ideokwargs.pop('capsize', None)
    ideokwargs.pop('orientation', None)

    fs = 5
    textkwargs = {
        "path_effects": [
            pe.withStroke(linewidth=fs / 3, foreground="white", alpha=0.8)
        ],
        "fontsize": fs,
        "va": "center_baseline",
        "ha": "center",
        "weight": "bold",
        "clip_on": True,
    } | textkwargs

    if names:
        for n, s, e in zip(band_name, band_start, band_end):
            midp = (e - s) // 2 + s
            ax.text(midp, 0, n, **textkwargs)

    figW, figH = ax.get_figure().get_size_inches()
    _, _, w, h = ax.get_position().bounds
    disp_ratio = (figW * w) / (figH * h)
    ratio = ax.get_data_ratio() * disp_ratio

    centro_pos = np.where(band_type == "acen")[0]
    round_patch = [
        barpatches[0],
        barpatches[band_num - 1],
    ]

    angled_patch = []

    if style == "angled":
        if orientation == "vertical":
            warnings.warn(
                "Angled look is not supported on vertical plots right now")
        angled_patch += [barpatches[centro_pos[0]], ax.patches[centro_pos[1]]]
    else:
        round_patch += [
            barpatches[centro_pos[1]],
            barpatches[centro_pos[0]],
        ]

    for i, patch in enumerate(round_patch):
        bb = patch.get_bbox()
        patch.remove()

        if orientation == "vertical":
            Box = SideRound(
                yround=min([((round_pct*0.9)*ratio) / bb.height, 2]),
                xround=0.9,
                corners=cornerref[i % 2],
            )
        else:
            Box = SideRound(
                xround=min([((round_pct*0.9)/ratio) / bb.width, 2]),
                yround=0.9,
                corners=cornerref[i % 2],
            )

        fp = FancyBboxPatch(
            (bb.xmin, bb.ymin),
            abs(bb.width),
            abs(bb.height),
            boxstyle=Box,
            ec=patch.get_edgecolor(),
            fc=patch.get_facecolor(),
            linewidth=patch.get_linewidth(),
            hatch=patch.get_hatch(),
            label=patch.get_label(),
        )

        fp._internal_update(ideokwargs)
        ax.add_patch(fp)

    for i, patch in enumerate(angled_patch):
        bb = patch.get_bbox()
        patch.remove()

        Box = [AngledBox(pad=0), ReverseFancyPatch(AngledBox(pad=0))][i % 2]

        fp = FancyBboxPatch(
            (bb.xmin, bb.ymin),
            abs(bb.width),
            abs(bb.height),
            boxstyle=Box,
            mutation_aspect=ratio,
            joinstyle="bevel",
            ec=patch.get_edgecolor(),
            fc=patch.get_facecolor(),
            linewidth=patch.get_linewidth(),
            hatch=patch.get_hatch(),
            label=patch.get_label(),
        )
        fp._internal_update(ideokwargs)
        ax.add_patch(fp)

    return ax


def ideogramv(
    chrom: str,
    bands=None,
    ax: plt.Axes = None,
    round_pct: float = 0.55,
    color: dict = BANDCOL,
    names=False,
    textkwargs={},
    **ideokwargs
):

    return ideogramh(
        chrom,
        bands=bands,
        ax=ax,
        round_pct=round_pct,
        color=color,
        style="round",
        names=names,
        textkwargs=textkwargs,
        orientation="vertical",
        ** ideokwargs
    )


def full_ideogram(bands=None, fig=None, orientation="h"):
    """
    This generates a full ideogram for the specified bands or for hg38,
    if bands are 'None'. I made this ages ago and testing it still work.
    Only 24 chromosomes are supported right now.

    Parameters
    ----------
    bands: filepath, dict or None, optional
        Bands dict. Need to be loaded already. The default is None.
    fig: optional
        A matplotlib figure object. If none is given, a new one will be made.
    orientation: str, optional
        Orientation. 'vertical' | 'v' or 'horizontal' | 'h' are supported.
        The default is "h".

    Raises
    ------
    ValueError
        If wrong orientation is given.

    Returns
    -------
    fig: TYPE
        The matplotlib figure that was made in this function.

    """

    # TODO: Load bands file and also support non-human bands

    if bands is None:
        with pkg_resources.as_file(
            pkg_resources.files(annotation) / "cytoBands.pickle.xz"
        ) as f:
            bands = pickle.load(lzma.open(f))
        chrnum = 24
        chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"] + ["chrY"]
    else:
        if isinstance(bands, dict):
            pass
        else:
            try:
                bands = load_cytobands(bands)
            except:
                ValueError("Bands must be a dict, file or None")

        chrnum = len(bands)
        chroms = list(bands.keys())

    if orientation not in {"h", "horizontal", "v", "vertical"}:
        raise ValueError(
            "Orientation needs to be either horizontal or vertical"
        )
    h = False
    if orientation in {"h", "horizontal"}:
        h = True

    if fig == None:
        fig = plt.figure(figsize=(20, chrnum) if h else (chrnum, 20))

    axes = []
    for i, chrom in enumerate(chroms):
        kwargs = (
            {"sharey" if orientation == "v" else "sharex": axes[0]}
            if i != 0
            else {}
        )

        if orientation == "h":
            ax = fig.add_subplot(chrnum, 1, i + 1, **kwargs)
        else:
            ax = fig.add_subplot(1, chrnum, i + 1, **kwargs)

        axes.append(ax)

        ax.xaxis.set_tick_params(
            which="both", length=0, **({"labelleft": False} if h else {})
        )
        ax.yaxis.set_tick_params(
            which="both", length=0, **({"labelleft": False} if not h else {})
        )
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["top"].set_visible(False)

        match orientation:
            case "h":
                ideogramh(chrom, bands, ax)
            case "v":
                ideogramv(chrom, bands, ax)

    return fig


def zoom(
    ax1,
    ax2,
    /,
    limits=None,
    zoomtype="curve",
    bend = 0.1,
):
    """
    Plot the ideogram of a chromosome

    Parameters
    ----------
    ax1 : plt.Axes
        Mandatory. The first axis to connect from.
    ax2 : plt.Axes
        Mandatory. The second axis to connect to.
    limits : plt.Axes
        x limits for zoom
    zoomtype : plt.Axes
        allowed are: "curve", "filledcurve". Anything else will make just a line
    bend : plt.Axes
        bend of the bezier curve. If line this is ingnored.
    """

    # TODO: allow connection from Y
    # TODO: add rectangle and connection properties
    # TODO: split fill from type!
    if (fig := ax1.get_figure()) != ax2.get_figure():
        raise ValueError("Axes need to be from the same figure!")

    if limits == None:
        start, end = ax1.get_xlim()
    else:
        try:
            start, end = limits
        except:
            raise ValueError(
                "limits needs to be either a tuple of integers that define the \
                limits of the zoomed plot or 'None' to have them \
                automatically determined"
            )

    if ax2.get_position().y0 > ax1.get_position().y0:
        ax1, ax2 = ax2, ax1

    y1lo, y1hi = ax1.get_ylim()
    y2lo, y2hi = ax2.get_ylim()

    # if rectangle == "ax2" or rectangle == "both":
    #     ax2.add_patch(
    #         plt.Rectangle(
    #             (start, y2lo), end - start, y2hi - y2lo, ec="k", alpha=0.6
    #         )
    #     )

    # if rectangle == "ax1" or rectangle == "both":
    #     ax1.add_patch(
    #         plt.Rectangle(
    #             (start, y1lo), end - start, y1hi - y1lo, ec="k", alpha=0.6
    #         )
    #     )

    if ax1.get_xticklabels():
        y1 = ax1.transData.inverted().transform(
            ax1.get_xticklabels()[0].get_window_extent().bounds[:2]
        )[1]
        y1 = min(y1, y1lo)
    else:
        y1 = y1lo

    y2 = y2hi

    # abs for inverted axes
    # get the absolute eight of the subplot. The whole y axis would be ~that in
    # percentage of the whole figure. Divide by <bend> to get the ratio of that
    # from 10% and divide the y-axis height to get the yaxis that is 10% of the
    # whole figure
    ax1p = abs(y1hi - y1lo) / ((ax1.bbox.height / fig.bbox.height) / bend)
    ax2p = abs(y2hi - y2lo) / ((ax2.bbox.height / fig.bbox.height) / bend)

    FILL = False
    if zoomtype == "curve" or zoomtype == "filledcurve":
        POINTS = [
            (start, y2),
            (start, y2 + ax2p),
            (start, y1 - ax1p),
            (start, y1),
            (end, y1),
            (end, y1 - ax1p),
            (end, y2 + ax2p),
            (end, y2),
        ]
        CODES = [
            mpath.Path.MOVETO,
            mpath.Path.CURVE4,
            mpath.Path.CURVE4,
            mpath.Path.CURVE4,
            mpath.Path.MOVETO,
            mpath.Path.CURVE4,
            mpath.Path.CURVE4,
            mpath.Path.CURVE4,
        ]
        COORDS = [
            ax2.transData,
            ax2.transData,
            ax1.transData,
            ax1.transData,
            ax1.transData,
            ax1.transData,
            ax2.transData,
            ax2.transData,
        ]

        if zoomtype == "filledcurve":
            CODES = [
                mpath.Path.MOVETO,
                mpath.Path.CURVE4,
                mpath.Path.CURVE4,
                mpath.Path.CURVE4,
                mpath.Path.LINETO,
                mpath.Path.CURVE4,
                mpath.Path.CURVE4,
                mpath.Path.CURVE4,
                mpath.Path.LINETO,
            ]
            COORDS.append(COORDS[0])
            POINTS.append(POINTS[0])
            FILL = True
    else:
        POINTS = [(start, y2), (start, y1), (end, y1), (end, y2)]
        CODES = [
            mpath.Path.MOVETO,
            mpath.Path.LINETO,
            mpath.Path.MOVETO,
            mpath.Path.LINETO,
        ]
        COORDS = [
            ax2.transData,
            ax1.transData,
            ax1.transData,
            ax2.transData,
        ]

    con = Chromosomeconnection(
        points=POINTS,
        path_codes=CODES,
        coords=COORDS,
        arrowstyle="-",
        shrinkB=5,
        facecolor=(0.1, 0.2, 0.5, 0.3),
        edgecolor="k",
        linewidth=1,
        fillable=FILL,
    )

    fig.add_artist(con)
