import shutil

import pytest


@pytest.fixture(scope="session")
def save_path(tmpdir_factory):
    """Docstring for save_path."""
    dir = tmpdir_factory.mktemp("temp_data", numbered=False)
    path = str(dir)
    yield path + "/"
    shutil.rmtree(str(tmpdir_factory.getbasetemp()))
