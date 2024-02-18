#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:15:09 2023

@author: balthasar
"""

import matplotlib.patches as mpatches
import numpy as np
from matplotlib.path import Path
import matplotlib.transforms as transforms
import matplotlib.path as mpath
from matplotlib._enums import CapStyle


# from matplotlib.path import Path as pltPath

# This class is the connection arrow, but works for other shapes as well
# I did not find an easier way to change it than to copy the whole code and
# and change many things inbetween.
# This is a problem for the license


class Chromosomeconnection(mpatches.FancyArrowPatch):
    """
    This is an variation of fancy arrow and connectionpatch to connect two plots
    with two lines.
    """

    _edge_default = True

    def __str__(self):
        return "ConnectionPatch((%g, %g), (%g, %g))" % (
            self.xy1[0],
            self.xy1[1],
            self.xy2[0],
            self.xy2[1],
        )

    def __init__(
        self,
        points,
        coords,
        path_codes=None,
        fillable=False,
        **kwargs,
    ):

        # if coordsB is None:
        #    coordsB = coordsA
        # we'll draw ourself after the artist we annotate by default
        self.points = points
        self.coords = coords
        self.pcodes = path_codes
        self.fillable = fillable

        # if True, draw annotation only if self.xy is inside the axes
        self._annotation_clip = None

        super().__init__(posA=(0, 0), posB=(1, 1), **kwargs)

        # kwargs.setdefault("joinstyle", JoinStyle.round)
        kwargs.setdefault("capstyle", CapStyle.butt)

    def _get_path_in_displaycoord(self):
        """Return the mutated path in display coordinates."""
        # TODO dpi correction?
        # see here https://github.com/matplotlib/matplotlib/blob/v3.8.2/lib/matplotlib/patches.py#L4054
        # how this is handled by matplotlib
        pos = [c.transform(xy) for c, xy in zip(self.coords, self.points)]

        path = [mpath.Path(pos, self.pcodes)]
        fillable = [self.fillable]

        return path, fillable

    def draw(self, renderer):
        if renderer is not None:
            self._renderer = renderer
        if not self.get_visible() or not self._check_xy(renderer):
            return

        self._dpi_cor = renderer.points_to_pixels(1.0)
        path, fillable = self._get_path_in_displaycoord()

        if not np.iterable(fillable):
            path = [path]
            fillable = [fillable]

        affine = transforms.IdentityTransform()

        self._draw_paths_with_artist_properties(
            renderer,
            [
                (
                    p,
                    affine,
                    self._facecolor if f and self._facecolor[3] else None,
                )
                for p, f in zip(path, fillable)
            ],
        )

    def set_annotation_clip(self, b):
        self._annotation_clip = b
        self.stale = True

    def get_annotation_clip(self):
        return self._annotation_clip

    def _check_xy(self, renderer):
        # TODO: I forgot what this is
        return True


class ReverseFancyPatch:
    """This is a class that given any FancyBboxPatch it will horizontally flip
    the returned Path"""

    def __init__(self, initialized_patch):

        self.patch = initialized_patch

    def __call__(self, x0, y0, width, height, mutation_size, **kwargs):

        # TODO: can I hide the wrapping from error messages inside patch?
        p = self.patch(x0, y0, width, height, mutation_size, **kwargs)
        p.vertices[:, 0] = 2 * x0 + width - p.vertices[:, 0]
        return p

# TODO: add pad and handle it better than I do right now.
# The fundamental problem with pad is that it changes x and y coords of the
# corners, but is necessary for the looks of the different shapes.
# A more general approach to this problem would be nice.


class SideRound:
    """A box with rounded corners."""

    def __init__(self, xround=0, yround=0, corners=[]):
        """
        Init.

        Parameters
        ----------
        xround : float
            amount of horizontal rounding

        yround : float
            amount of vertical rounding

        corners : list
            list of corners that should be rounded. 0 is the bottom left and
            from there it is counted counter-clockwise.
        """

        self.xround = xround
        self.yround = yround
        self.corners = corners
        super().__init__()

    def __call__(self, x0, y0, width, height, mutation_size):

        xdr = self.xround
        ydr = self.yround

        x1 = x0 + width
        y1 = y0 + height

        drx = ((x1 - x0) / 2) * xdr
        dry = ((y1 - y0) / 2) * ydr

        ROUND = [
            [(x0, y0 + dry), (x0, y0), (x0 + drx, y0)],
            [(x1 - drx, y0), (x1, y0), (x1, y0 + dry)],
            [(x1, y1 - dry), (x1, y1), (x1 - drx, y1)],
            [(x0 + drx, y1), (x0, y1), (x0, y1 - dry)],
        ]

        EDGY = [
            [(x0, y0)],
            [(x1, y0)],
            [(x1, y1)],
            [(x0, y1)],
        ]

        cp = []
        com = [mpath.Path.MOVETO]

        for i in range(4):
            if i in self.corners:
                cp += ROUND[i]
                com += [
                    mpath.Path.CURVE3,
                    mpath.Path.CURVE3,
                    mpath.Path.LINETO,
                ]
            elif i == 0 or cp[-1] != EDGY[i]:
                cp += EDGY[i]
                com += [mpath.Path.LINETO]

        if cp[-1] == cp[0] and com[-2] == mpath.Path.LINETO:
            com.pop()
            com[-1] = mpath.Path.CLOSEPOLY
        else:
            cp.append(cp[0])
            com[-1] = mpath.Path.CLOSEPOLY

        path = mpath.Path(cp, com)

        return path


class AngledBox:
    """An angled box."""

    def __init__(self, pad=0.3, angle=45):
        """
        The arguments must be floats and have default values.

        Parameters
        ----------
        pad : float
            amount of padding
        angle : float
            angle of the arrow side of the box. default is 45
        """
        self.pad = pad
        self.angle = np.deg2rad(angle)
        super().__init__()

    def __call__(self, x0, y0, width, height, mutation_size):
        """
        Given the location and size of the box, return the path of the box
        around it.

        Rotation is automatically taken care of.

        Parameters
        ----------
        x0, y0, width, height : float
            Box location and size.
        mutation_size : float
            Reference scale for the mutation, typically the text font size.
        """
        # padding
        pad = mutation_size * self.pad
        # width and height with padding added
        width = width + 2.0 * pad
        height = height + 2.0 * pad

        angle = self.angle
        # height is twice the gegenkathete and we need the ankathete
        margin = height / (np.tan(angle) * 2)

        if margin > width:
            margin = width
        # boundary of the padded box
        x0, y0 = x0 - pad, y0 - pad
        x1, y1 = x0 + width, y0 + height
        # return the new path

        old = (None, None)
        path = [
            (old := point)
            for point in [
                (x0, y0),
                (x0, y1),
                (x1 - margin, y1),
                (x1 + pad, (y0 + y1) / 2.0),
                (x1 - margin, y0),
                (x0, y0),
            ]
            if point != old
        ]

        return Path(path, closed=True)


class ThinArrow:
    """An flat arrow"""

    def __init__(self, pad=0.3, angle=45):
        """
        The arguments must be floats and have default values.

        Parameters
        ----------
        pad : float
            amount of padding
        angle : float
            angle of the arrow in degree. default is 45
        """
        self.pad = pad
        self.angle = np.deg2rad(angle)
        super().__init__()

    def __call__(self, x0, y0, width, height, mutation_size):
        """
        Given the location and size of the box, return the path of the box
        around it.

        Rotation is automatically taken care of.

        Parameters
        ----------
        x0, y0, width, height : float
            Box location and size.
        mutation_size : float
            Reference scale for the mutation, typically the text font size.
        """
        # padding
        pad = mutation_size * self.pad
        # width and height with padding added
        width = width + 2.0 * pad
        height = height + 2.0 * pad

        angle = self.angle
        # height is twice the gegenkathete and we need the ankathete
        margin = height / (np.tan(angle) * 2)

        if margin > width:
            margin = width
        # boundary of the padded box
        x0, y0 = x0 - pad, y0 - pad
        x1, y1 = x0 + width, y0 + height
        # return the new path

        return Path(
            [
                (x0, y0),
                (x0, y1),
                (x0, y0 + height / 2),
                (x1, y0 + height / 2),
                (x1 - margin, y1),
                (x1 + pad, (y0 + y1) / 2.0),
                (x1 - margin, y0),
                # (x0, y0),
            ],
            codes=[
                mpath.Path.MOVETO,
                mpath.Path.LINETO,
                mpath.Path.MOVETO,
                mpath.Path.LINETO,
                mpath.Path.MOVETO,
                mpath.Path.LINETO,
                mpath.Path.LINETO,
            ],
        )


class Arrowed:
    """An arrowed line."""

    def __init__(self, pad=0.3, angle=45, arrowscale=0.5):
        """
        The arguments must be floats and have default values.

        Parameters
        ----------
        pad : float
            amount of padding
        angle : float
            angle of the arrow in degree. default is 45
        arrowscale: float
            I how high a single wing of the arrow should be as fraction of height.
            0.5 (the default) means the arrow spans the whole height.
        """
        self.pad = pad
        self.angle = np.deg2rad(angle)
        self.arrowscale = arrowscale
        super().__init__()

    def __call__(self, x0, y0, width, height, mutation_size):
        """
        Given the location and size of the box, return the path of the box
        around it.

        Rotation is automatically taken care of.

        Parameters
        ----------
        x0, y0, width, height : float
            Box location and size.
        mutation_size : float
            Reference scale for the mutation, typically the text font size.
        """
        # padding
        pad = mutation_size * self.pad
        # width and height with padding added
        width = width + 2.0 * pad
        height = height + 2.0 * pad

        angle = self.angle
        # height is twice the gegenkathete and we need the ankathete
        aheight = (height - height * self.arrowscale) / 2

        margin = aheight / (np.tan(angle) * 2)

        if margin > width:
            margin = width
        # boundary of the padded box
        x0, y0 = x0 - pad, y0 - pad
        x1, y1 = x0 + width, y0 + height
        # return the new path

        mid = np.linspace(x0, x1, max(2, round(width / (2 * height))))
        top = mid - margin

        xs = [val for pair in zip(top, mid, top) for val in pair]
        ys = [y0 + aheight, y0 + height / 2, y1 - aheight] * len(mid)

        return Path(
            [
                (x1, y0 + height / 2),
                (x0, y0 + height / 2),
                *zip(xs, ys),
            ],
            codes=[
                mpath.Path.MOVETO,
                mpath.Path.LINETO,
                *(
                    (mpath.Path.MOVETO, mpath.Path.LINETO, mpath.Path.LINETO)
                    * len(mid)
                ),
            ],
        )
