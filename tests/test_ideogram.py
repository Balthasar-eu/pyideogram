import pytest
import matplotlib
matplotlib.use("Agg")  # Use a non-interactive backend for testing
from matplotlib import pyplot as plt

import sys
import os
# Ensure import works if running from repo root
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from pyideogram.ideogram import ideogramh, ideogramv

@pytest.mark.parametrize("func,orientation", [
    (ideogramh, "horizontal"),
    (ideogramv, "vertical"),
])
def test_ideogram_functions(func, orientation):
    chrom = "chr6"
    fig, ax = plt.subplots(figsize=(6, 1 if orientation == "horizontal" else 6))
    result_ax = func(chrom, ax=ax)
    assert result_ax is ax
    assert hasattr(ax, "patches")
    assert len(ax.patches) > 0

def test_ideogramh_default():
    # Test ideogramh with minimal args
    chrom = "chr1"
    fig, ax = plt.subplots(figsize=(6, 1))
    out = ideogramh(chrom, ax=ax)
    assert out is ax
    assert hasattr(ax, "patches")
    assert len(ax.patches) > 0

def test_ideogramv_default():
    # Test ideogramv with minimal args
    chrom = "chr1"
    fig, ax = plt.subplots(figsize=(1, 6))
    out = ideogramv(chrom, ax=ax)
    assert out is ax
    assert hasattr(ax, "patches")
    assert len(ax.patches) > 0
