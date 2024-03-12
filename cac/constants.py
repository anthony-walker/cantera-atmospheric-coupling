import os

VERSION = '0.0.1'
SRC_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(SRC_DIR, 'data')
TEST_DIR = os.path.join(SRC_DIR, "tests")
TEST_DATA_DIR = os.path.join(TEST_DIR, "data")
COLORS = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']
MARKERS = ['o', 's', '^', 'v', '>', '<', 'x', '+', '*', 'D', 'd', 'h', 'H', 'p', '|', '_']
HATCHES = [
    '/', '\\', '-', 'x', '//', '|', '+', '\\', '*', 'o', 'O', '.', '-',
    '\\|', '//|', '-', 'v', '^', '<', '>', '1', '2', '3', '4'
]
