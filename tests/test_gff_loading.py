import os
import pytest

# Adjust path for local import if needed
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from pyideogram.dataloader import load_gff_tree, GeneTree

def get_test_file(filename):
    return os.path.join(os.path.dirname(__file__), "data", filename)

def test_load_gff_tree_dict():
    gff_path = get_test_file('test.gff')
    trees = load_gff_tree(gff_path)
    assert isinstance(trees, dict)
    assert trees  # not empty
    for chrom, tree in trees.items():
        assert isinstance(tree, GeneTree)

def test_load_gff_tree_single_chrom():
    gff_path = get_test_file('test.gff')
    # Assume one of the chromosomes is 'chr1' in your test.gff
    tree = load_gff_tree(gff_path, chromosome="chr1")
    assert isinstance(tree, GeneTree)
