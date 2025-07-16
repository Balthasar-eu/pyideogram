#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 12:26:18 2023

@author: balthasar
"""


import warnings
import bisect

###############################################################################
# %% ################################ classes ################################
###############################################################################

# The design of these classes might seem strange, but it is to save them efficiently
# as pickles and still have dict like access to their attributes.


class Exon:
    __slots__ = "start", "end"

    def __repr__(self):
        return f"{self.__class__.__name__}({self.start}, {self.end})"

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __lt__(self, other):
        if self.start > self.end:
            if other.start < other.end:
                raise ValueError("Cannot mix forward and reverse Exons")
            return -self.end < -other.end
        else:
            return self.start < other.start


class Transcript:
    __slots__ = "name", "exons", "__dict__"

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.name}', {self.exons})"

    def __init__(self, name, exons, **kwargs):
        self.name = name
        self.exons = sorted(exons)
        self.__dict__ = kwargs

    def __getitem__(self, key):
        match key:
            case "name":
                return self.name
            case "exons":
                return self.exons
            case "start":
                if key in self.__dict__:
                    return self.__dict__[key]
                else:
                    return self.exons[0].start
            case "end":
                if key in self.__dict__:
                    return self.__dict__[key]
                else:
                    return self.exons[-1].end
            case _:
                if key in self.__dict__:
                    return self.__dict__[key]
                else:
                    raise ValueError(
                        f"The transcript has no attribute '{key}'!"
                    )

    def __setitem__(self, key, value):
        match key:
            case "name":
                self.name = value
            case "exons":
                self.exons = value
            case _:
                self.__dict__[key] = value


class Gene:
    """
    This is a class for gene info. Has slots for name, start and end and transcripts.
    """

    __slots__ = "name", "start", "end", "transcripts", "__dict__"

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.name}', {self.start}, {self.end}, {self.transcripts})"

    def __init__(self, name, start, end, transcripts, **kwargs):
        self.name = name
        self.transcripts = transcripts
        self.start = start
        self.end = end
        self.__dict__ = kwargs

    def __getitem__(self, key):
        match key:
            case "name":
                return self.name
            case "transcripts":
                return self.transcripts
            case _:
                if key in self.__dict__:
                    return self.__dict__[key]
                else:
                    raise ValueError(f"The gene has no attribute '{key}'!")

    def __lt__(self, other):
        return self.start < other.start

    def __setitem__(self, key, value):
        match key:
            case "name":
                self.name = value
            case "transcripts":
                self.transcripts = value
            case _:
                self.__dict__[key] = value


class GeneTree:
    """
    Not a tree, but it sorts a gene list and performs binary search. This works
    well in the specific usecase of genes, because genes are generally much
    smaller than the search space and non overlapping. I found this can perform two orders of magnitude
    better than intervaltree for specific queries.
    """

    __slots__ = "genes", "largest"

    def __init__(self, genes):
        if not isinstance(genes, list):
            raise ValueError("'genes' must be a list")

        self.genes = sorted(genes)
        self.largest = max([g.end - g.start for g in genes])

    def __getitem__(self, index):
        start, stop = index.start, index.stop

        if start is None:
            start = 0
            if stop is None:
                return []
        if stop is None:
            stop = self.genes[-1].end + 1

        mins = max(start - self.largest, 0)
        i = bisect.bisect_left(self.genes, Gene("dummy", mins, mins + 1, []))

        glist = []
        for g in self.genes[i:]:
            if g.start > stop:
                break
            if g.end > start:
                glist.append(g)

        return glist


def load_cytobands(file):
    banddict = {}
    with open(file) as f:
        for i, line in enumerate(f):
            try:
                chrom, start, end, name, btype = line.strip().split()
            except:
                warnings.warn(f"Error in line {i}. Skipping.")
            else:
                blists = banddict.setdefault(chrom, ([], [], []))
                blists[0].append(int(end))
                blists[1].append(name)
                blists[2].append(btype)
    return banddict


def load_gff_tree(
    file,
    chromosome=None,
    name_field="gene_name",
    tname_field="transcript_name",
):
    """

    This loads a gff into a tree like structure to parse genes quickly by pos.
    It is a bit untested and WIP.

    Parameters
    ----------
    file : str or path
        The path to the gff
    chromosome : str or None
        The chromosome to load, or leave empty for everything

    Returns
    -------
    TYPE
        A dictionary of chromosomes or a 'Tree' of all genes on a specifc chromosome.

    """

    # TODO clean up the whole function and test that everything works well with
    # the new classes!

    # fmt: off

    names = [
        "chrom", "source", "feature", "start", "end", "score", "strand", "phase",
    ]

    types = [ str, str, str, int, int, str, str, str ]

    # fmt: on

    # TODO: split everything by chromosome

    transcripts = {}
    exons = {}
    genes = {}
    with open(file) as f:
        for line in f:
            stripline = line.strip().split("\t")
            if stripline[0][0] == "#":
                continue

            cand = {
                n: t(s) for n, s, t in zip(names, stripline[:-1], types)
            } | dict(
                [i.strip('"') for i in item.strip().split("=", 1)]
                for item in stripline[-1].split(";")
                if item.strip()
            )

            # only include transcripts and exons
            match cand["feature"]:
                case "gene":
                    genes.setdefault(cand["chrom"], []).append(cand)
                case "transcript":
                    transcripts.setdefault(cand["chrom"], []).append(cand)
                case "exon":
                    exons.setdefault(cand["chrom"], []).append(cand)
    gdicts = {}
    errorcount = 0
    for chrom in genes:
        tdict = {}
        if chromosome is not None and chromosome != chrom:
            continue
        gdict = {}
        for g in genes[chrom]:
            if g["ID"] in gdict:
                print("Duplicate gene ID", g["ID"])
            else:
                try:
                    gdict[g["ID"]] = (
                        g[name_field],
                        g["start"],
                        g["end"],
                        g["strand"],
                        {},
                    )
                except KeyError:
                    errorcount += 1
                    continue

        for t in transcripts[chrom]:
            try:
                parent = t["Parent"]
                tdict[t["ID"]] = t["Parent"]
                gdict[parent][4][t["ID"]] = {
                    "name": t[tname_field],
                    "exons": [],
                }
            except KeyError:
                errorcount += 1
                continue

        for e in exons[chrom]:
            try:
                gparent = tdict[e["Parent"]]
                tparent = e["Parent"]
                gdict[gparent][4][tparent]["exons"].append(
                    (
                        e["start"],
                        e["end"],
                    )
                )
            except KeyError:
                continue
        gdicts[chrom] = gdict

    if errorcount:
        warnings.warn(
            f"""Encountered {errorcount} key errors. Make sure the name field
            in the gff is {name_field}. If it is not you can change it by
            giving the correct name as function argument."""
        )

    for chrom, gdict in gdicts.items():
        genes = []

        for g, v in gdict.items():
            if v[3] == "-":
                new = Gene(
                    v[0],
                    v[1],
                    v[2],
                    [
                        Transcript(
                            tattr["name"],
                            [Exon(tte, tts) for tte, tts in tattr["exons"]],
                        )
                        for t, tattr in v[4].items()
                    ],
                )
            else:
                new = Gene(
                    v[0],
                    v[1],
                    v[2],
                    [
                        Transcript(
                            tattr["name"],
                            [Exon(tts, tte) for tts, tte in tattr["exons"]],
                        )
                        for t, tattr in v[4].items()
                    ],
                )

            genes.append(new)

        gdicts[chrom] = GeneTree(genes)

    if chromosome is None:
        return gdicts
    else:
        return gdicts[chromosome]


try:
    import pandas as pd
except ImportError:
    pass
else:

    def load_gff_pd(file, limits=None, ordered=False):
        """

        THIS IS CURRENTLY NOT USED, BUT SHOULD WORK.
        This approach is somewhat manual, but faster than splitting and
        converting columns after loading the file.

        Parameters
        ----------
        file : TYPE
            DESCRIPTION.
        limits : TYPE, optional
            If coordinates are outside this range, we don't load them.
            The default is None.
        ordered : TYPE, optional
            Is the gff ordered? Only relevant if limits are given, because
            then we can skip parsing the whole file. The default is False.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """

        # ordered will end the loading, when a value higher than the upper limit is
        # encountered

        names = [
            "chrom",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "phase",
        ]

        types = [
            str,
            str,
            str,
            int,
            int,
            str,
            str,
            str,
        ]

        ordered = ordered and limits is None

        if limits is None:
            unlimited = True
            start = end = None
        else:
            unlimited = False
            start, end = limits
        dfcols = []
        with open(file) as f:
            for line in f:
                stripline = line.strip().split("\t")
                if stripline[0][0] == "#":
                    continue
                if unlimited or (
                    int(stripline[4]) > start and int(stripline[3]) < end
                ):
                    dfcols.append(
                        {
                            n: t(s)
                            for n, s, t in zip(names, stripline[:-1], types)
                        }
                        | dict(
                            [i.strip('"') for i in item.strip().split("=", 1)]
                            for item in stripline[-1].split(";")
                            if item.strip()
                        )
                    )
                if ordered and int(stripline[3]) > end:
                    break

        return pd.DataFrame(dfcols)
