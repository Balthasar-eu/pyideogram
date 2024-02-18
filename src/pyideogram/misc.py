#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 18:20:10 2023

@author: balthasar
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

###############################################################################
#%% ################################ functions ################################
###############################################################################


def set_genome_xticks(ax=None, start=None, end=None, steps=10, prune=None):

    if ax == None:
        ax = plt.gca()

    tmp_start, tmp_end = ax.get_xlim()

    if start == None:
        start = tmp_start
    if end == None:
        end = tmp_end

    ax.xaxis.set_major_formatter(
        matplotlib.ticker.EngFormatter(unit="b", sep="")
    )
    # ax.xaxis.set_major_locator(
    #     matplotlib.ticker.MultipleLocator(
    #         base=10
    #         ** (
    #             int(np.log10(end - start))
    #             - (1 if int(str(end - start)[:2]) < 15 else 0)
    #             # if first digit is 1, only one tick would be generated
    #             # this makes it 10
    #         )
    #     )
    # )

    ax.xaxis.set_major_locator(
        matplotlib.ticker.MaxNLocator(
            steps=[steps], min_n_ticks=0, prune=prune
        )
    )

    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(10))
    # labels = ax.axes.get_xticklabels()
    # if len(labels) > 5:
    #    plt.setp(plot.axes.get_xticklabels(), visible=False)
    #    plt.setp(plot.axes.get_xticklabels()[::5], visible=True)


def set_genome_yticks(ax=None, start=None, end=None):

    if ax == None:
        ax = plt.gca()

    tmp_start, tmp_end = ax.get_ylim()

    if start == None:
        start = tmp_start
    if end == None:
        end = tmp_end

    ax.yaxis.set_major_formatter(
        matplotlib.ticker.EngFormatter(unit="b", sep="")
    )
    ax.yaxis.set_major_locator(
        matplotlib.ticker.MultipleLocator(
            base=10 ** (int(np.log10(end - start)))
        )
    )
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(10))
