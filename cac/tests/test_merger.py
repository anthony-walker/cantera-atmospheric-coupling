import os
import shutil
import pytest
import unittest
from cac.merger import *
from cac.constants import TEST_DATA_DIR

class TestMerger(unittest.TestCase):

    def setUp(self) -> None:
        pass

    # @pytest.mark.diagnose
    def test_mechanism_merge(self):
        fmech = os.path.join(TEST_DATA_DIR, "fuel.yaml")
        amech = os.path.join(TEST_DATA_DIR, "atms.yaml")
        merge_mechanisms(fmech, amech)

    @pytest.mark.diagnose
    def test_mechanism_update(self):
        update_file = os.path.join(TEST_DATA_DIR, "update.yaml")
        copied_file = os.path.join(TEST_DATA_DIR, "updated-update.yaml")
        update_mechanism_with_smiles_and_inchi(update_file)
        assert os.path.exists(copied_file)
        # os.remove(copied_file)
