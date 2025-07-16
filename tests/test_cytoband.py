import os
import pytest
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

# Adjust path for local import if needed
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from pyideogram.ideogram import load_cytobands, ideogramh

def get_test_file(filename):
    return os.path.join(os.path.dirname(__file__), "data", filename)

def test_load_cytobands_from_file():
    cytopath = get_test_file('cytoBand.tsv')
    bands = load_cytobands(cytopath)
    assert isinstance(bands, dict)
    assert len(bands) > 0
    assert 'chr1' in bands  # loose check for dict or tuple

def test_ideogramh_with_custom_bands():
    cytopath = get_test_file('cytoBand.tsv')
    bands = load_cytobands(cytopath)
    fig, ax = plt.subplots(figsize=(6, 1))
    out = ideogramh("chr1", bands=bands, ax=ax)
    assert out is ax
    assert len(ax.patches) > 0

def test_ideogramh_with_custom_bands_from_file():
    cytopath = get_test_file('cytoBand.tsv')
    fig, ax = plt.subplots(figsize=(6, 1))
    out = ideogramh("chr1", bands=cytopath, ax=ax)
    assert out is ax
    assert len(ax.patches) > 0
