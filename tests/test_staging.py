"""Cache staging: a cache path is either absent or complete."""

import pytest

from shachen.io.staging import staged_download


def test_cache_appears_only_when_the_block_finishes(tmp_path):
    """A reader that looks mid-download must not find the cache.

    This is the parallel case: several scenes of one event ask for the same
    MERRA-2 day file at once, and the writer's partial file used to be
    visible under the cache's own name.
    """
    cache = tmp_path / "merra2_ts_20210315.nc"
    with staged_download(cache) as (scratch, staged):
        staged.write_bytes(b"half")
        assert not cache.exists()
        assert scratch.is_dir()
        staged.write_bytes(b"whole")
    assert cache.read_bytes() == b"whole"


def test_failure_leaves_neither_cache_nor_scratch(tmp_path):
    """An interrupted download must not publish, and must not accumulate."""
    cache = tmp_path / "merra2_met_20210315.nc"
    with pytest.raises(RuntimeError):
        with staged_download(cache) as (_, staged):
            staged.write_bytes(b"partial")
            raise RuntimeError("granule search failed")
    assert not cache.exists()
    assert list(tmp_path.iterdir()) == []


def test_concurrent_stages_do_not_share_a_scratch_directory(tmp_path):
    """Two processes staging the same cache each get their own scratch.

    ``earthaccess.download`` names the granule after the remote object, so
    sharing one directory meant one caller deleting or renaming the other's
    in-flight file.
    """
    cache = tmp_path / "CAM5K30EM_202103.nc"
    with staged_download(cache) as (first_dir, first):
        with staged_download(cache) as (second_dir, second):
            assert first_dir != second_dir
            second.write_bytes(b"from the second")
        # The loser of the race publishes first; the winner overwrites it
        # with an equally complete file, which is what makes it harmless.
        assert cache.read_bytes() == b"from the second"
        first.write_bytes(b"from the first")
    assert cache.read_bytes() == b"from the first"


def test_creates_the_cache_directory(tmp_path):
    """The ancillary directories are created on first use, as before."""
    cache = tmp_path / "merra2" / "merra2_ts_20210315.nc"
    with staged_download(cache) as (_, staged):
        staged.write_bytes(b"x")
    assert cache.is_file()
