import shutil
from distutils.dir_util import copy_tree

import pytest


@pytest.fixture(scope="session")
def save_path(tmpdir_factory):
    """Docstring for save_path."""
    dir = tmpdir_factory.mktemp("temp_data", numbered=False)
    path = str(dir)
    copy_tree("tests/data", path)
    yield path + "/"
    shutil.rmtree(str(tmpdir_factory.getbasetemp()))
