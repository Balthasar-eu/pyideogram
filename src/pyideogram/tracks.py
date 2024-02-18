#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 18:20:10 2023

@author: balthasar
"""
import numpy as np
import matplotlib.pyplot as plt
from .dataloader import load_gff_tree, Exon, Transcript, Gene, GeneTree
from .matplotlib_extension import (
    ReverseFancyPatch,
    AngledBox,
    ThinArrow,
    Arrowed,
)
from matplotlib.patches import FancyBboxPatch
import matplotlib.patheffects as pe
from matplotlib.patches import BoxStyle

from pathlib import Path
import lzma
import pickle

from . import annotation
import importlib.resources as pkg_resources


_ANNOTATION_CACHE_ = {}
_ANNOTATION_BUILTIN_CACHE_ = {}


###############################################################################
# %% ################################ functions ################################
###############################################################################


# def bed(bed, chrom, ax=None, textcol=None, valuecol=None, ptype="plot"):
#     """

#     Not used right now

#     Parameters
#     ----------
#     bed : TYPE
#         DESCRIPTION.
#     chrom : TYPE
#         DESCRIPTION.
#     ax : TYPE, optional
#         DESCRIPTION. The default is None.
#     textcol : TYPE, optional
#         DESCRIPTION. The default is None.
#     valuecol : TYPE, optional
#         DESCRIPTION. The default is None.
#     ptype : TYPE, optional
#         DESCRIPTION. The default is "plot".

#     Raises
#     ------
#     ValueError
#         DESCRIPTION.

#     Returns
#     -------
#     ax : TYPE
#         DESCRIPTION.

#     """
#     if ax == None:
#         ax = plt.gca()

#     if type(bed) == str:
#         bed = pd.read_csv(bed, sep="\t", header=None, index_col=None)

#     if len(ax.get_shared_x_axes().get_siblings(ax)) > 1:
#         ax.xaxis.set_tick_params(which="both", length=0, labelbottom=False)

#     bed = bed.loc[bed[0] == chrom]
#     # needed for text and scatter
#     bed["midp"] = bed.iloc[:, 1] + (bed.iloc[:, 2] - bed.iloc[:, 1]) / 2

#     if valuecol == None:
#         ycol = [0] * len(bed)
#         ax.yaxis.set_tick_params(which="both", length=0, labelbottom=False)
#     else:
#         ycol = bed.iloc[:, valuecol]

#     if textcol:
#         if not ax.get_autoscale_on():
#             start, end = ax.get_xlim()
#             print(start, end)
#         inlimits = (bed["midp"] > start) & (bed["midp"] < end)
#         textcol = bed.iloc[:, textcol]
#         for m, t in zip(bed["midp"].loc[inlimits], textcol.loc[inlimits]):
#             ax.text(m, 0, t, ha="center", va="bottom")

#     if ptype == "plot":
#         for start, end, y in zip(bed.iloc[:, 1], bed.iloc[:, 2], ycol):
#             ax.plot([start, end], [0, 0])
#     elif ptype == "scatter":
#         ax.scatter(bed["midp"], ycol)
#     else:
#         raise ValueError(f"{ptype} is not an appropriate plot type")

#     return ax


###############################################################################
# %% ########################### plotting functions ###########################
###############################################################################


def genetrack(
    region,
    /,
    geneinfo=None,
    ax=None,
    lanenum=0,
    textlane=None,
    print_text=True,
    text_fit_height=True,
    text_clip_width=True,
    transcriptstyle=None,
    exonstyle=None,
    trackkwargs={},
    textkwargs={},
):
    """
    Plot gene annotation from a gff or gtf file

    Parameters
    ----------
    region :
        Mandatory. The genomic location to plot. Uses standard genomic
        coordinate format. <contig>:<start>-<end>
    geneinfo : file, str, pathlib.Path, intervaltree, optional
        A path to an annotation file in gff or gtf format or a
        intervaltree that was created from the load_gff/gtf function from this
        library. If no value is provided, the MANE annotations for hg38, which
        are included in this package, are used.
    ax : matplotlib.Axes, optional
        An matplotlib axes object that will be used to plot the genes, if
        provided. If not provided plt.gca will be called to get one.
    lanenum: int, optional
        number of lanes to plot genes one. If set to 0 or a negative number
        the number of lanes will be automatically determined in a way that
        guarantees no collisions will occur. If it is given, but too small,
        collisions might not be avoided. If a large number of genes is plotted,
        giving the number of lanes can speed up the plotting.
    text_fit_height: bool,optional
        Should the text inside genes be fitted to the height of the lane.
    text_clip_width: bool, optional
        should text be clipped if it is outside the bounding box of the plot
    boxstyle: Boxstyle, optional
        Boxstyle for the appearance of the transcripts/exons
    trackkwargs: dict, optional
        arguments for the underlying matplotlib functions that create the
        elements that compose the track.
    textkwargs: dict, optional
        arguments for the underlying matplotlib functions that create the
        text.

    Raises
    ------
    ValueError
        If the region does not conform to the established standard format.

    Returns
    -------
    None.


    """

    # This is arbitary, but is only used for display purposes, so no need to
    # have it user settable
    LANEHEIGHT = 0.7

    global _ANNOTATION_CACHE_

    ###########################################################################
    ############################# input  handling #############################

    if ax == None:
        ax = plt.gca()

    if len(ax.get_shared_x_axes().get_siblings(ax)) > 1:
        ax.xaxis.set_tick_params(which="both", length=0, labelbottom=False)

    # TODO: test if at least the first input is valid
    if isinstance(region, list):
        genes = region
    else:
        try:
            chrom, stend = region.split(":")
            start, end = stend.split("-")
            start = int(start)
            end = int(end)
        except ValueError:
            raise ValueError(
                """The region needs to conform to the standard genomic
                coordinate format: <CONTIG>:<START>-<END>"""
            )

        ###########################################################################
        ############################# annot  handling #############################

        if geneinfo is None:
            if chrom in _ANNOTATION_BUILTIN_CACHE_:
                geneinfo = _ANNOTATION_BUILTIN_CACHE_[chrom]
            else:
                annot_file = (
                    pkg_resources.files(annotation)
                    / f"MANE38_{chrom}.pickle.xz"
                )

                with pkg_resources.as_file(annot_file) as f:
                    if f.is_file():
                        geneinfo = pickle.load(lzma.open(f))
                    else:
                        raise ValueError(
                            f"""No default annotation available for contig "{chrom}"
                            Make sure that the contig name in the region argument is
                            correct, and/or provide an annotation gff/gtf.
                            """
                        )

                _ANNOTATION_BUILTIN_CACHE_[chrom] = geneinfo

        elif isinstance(geneinfo, GeneTree):
            print("tree given")
        elif isinstance(geneinfo, str) or isinstance(geneinfo, Path):
            geneinfo = load_gff_tree(geneinfo, chrom)
        else:
            print(
                """If any input for geneinfo is given, it needs to be either an
                  Intervaltree or a Path like object."""
            )

        genes = sorted(geneinfo[start:end])

    if ax.get_autoscale_on():
        ax.set_xlim(start, end)

    ###########################################################################
    ############################# calculate lanes #############################

    lanescale = 1
    if lanenum < 1:
        lanes = [0]
        for gene in genes:
            for transcript in gene["transcripts"]:
                tstart, tend = sorted([transcript["start"], transcript["end"]])
                if tstart < min(lanes):
                    lanes.append(tend)
                else:
                    lanes[np.argmin(lanes)] = tend
        lanenum = len(lanes)

    if textlane or (textlane is None and end - start > 200000):
        lanescale = 2
        textlane = True

    ax.set_ylim(-1, lanenum * lanescale)
    ax.set_yticks([])  # [(lanenum - 1) / 2])
    ax.set_yticklabels([])

    ###########################################################################
    ############################## get  fontsize ##############################
    # The initial fontsize is irrelevant, but it needs to be set to calculate
    # the optimal one.

    if print_text:
        fs = 5
        if text_fit_height:
            t = ax.text(
                0,
                0,
                " ",
                fontsize=fs,
                va="center_baseline",
                ha="center",
                weight="bold",
            )

            # Scale to height
            bbox = t.get_window_extent().transformed(ax.transData.inverted())
            fs = fs * (LANEHEIGHT * 0.7) / bbox.height
            t.set_visible(False)

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

    ###########################################################################
    ############################### pixel ratio ###############################
    # kinda complicated, but matplotlibs fault. See linked stackoverflow

    # ugh https://stackoverflow.com/questions/41597177/get-aspect-ratio-of-axes
    figW, figH = ax.get_figure().get_size_inches()
    _, _, w, h = ax.get_position().bounds
    disp_ratio = (figW * w) / (figH * h)

    # this should make the arrow angle 45 degree in most cases
    # the arrow should be as far as high. because LANEHEIGHT is double the arrow
    # height, get the pixel ratio and then divide by additional 2
    ratio = ax.get_data_ratio() * disp_ratio

    ###########################################################################
    #################################  kwargs #################################
    # These are some defaults but can be overwritten by the users.

    trackkwargs = {
        "ec": "black",
        "joinstyle": "round",
        "capstyle": "round",
        "color": "k",
    } | trackkwargs

    ###########################################################################
    ################################ plotting  ################################

    lanes = [0] * lanenum
    for gene in genes:
        tnum = len(gene["transcripts"])
        for transcript in gene["transcripts"]:
            tstart, tend = sorted([transcript["start"], transcript["end"]])
            tlen = tend - tstart

            y = np.argmin(lanes)
            lanes[y] = max(tend, lanes[y])

            y = (y * lanescale) - LANEHEIGHT / 2

            if callable(transcriptstyle):
                style = transcriptstyle
            else:
                match transcriptstyle:
                    case "arrow":
                        style = ThinArrow(pad=0)
                    case "angled":
                        style = AngledBox(pad=0)
                    case "arrowed":
                        style = Arrowed(pad=0)
                    case _:
                        if end - start > 100000:
                            style = AngledBox(pad=0)
                        else:
                            style = Arrowed(pad=0)

                if transcript["start"] > transcript["end"]:
                    style = ReverseFancyPatch(style)

                ax.add_patch(
                    FancyBboxPatch(
                        (tstart, y),
                        tlen,
                        LANEHEIGHT,
                        boxstyle=style,
                        mutation_aspect=ratio,
                        **trackkwargs,
                    )
                )

            if exonstyle is not None:
                for exon in transcript["exons"]:
                    exonlen = exon.end - exon.start
                    ax.add_patch(
                        FancyBboxPatch(
                            (exon.start, y),
                            exonlen,
                            LANEHEIGHT,
                            boxstyle=BoxStyle.Square(pad=0),
                            mutation_aspect=ratio,
                            **trackkwargs,
                        )
                    )

            if print_text:
                if tnum == 1:
                    textv = gene["name"]
                else:
                    textv = transcript["name"]

                if textlane:
                    y += 1

                tt = ax.text(
                    tstart + tlen / 2, y + LANEHEIGHT / 2, textv, **textkwargs
                )
                if (
                    text_clip_width
                    and tt.get_window_extent()
                    .transformed(ax.transData.inverted())
                    .width
                    > tlen
                ):
                    tt.set_visible(False)
