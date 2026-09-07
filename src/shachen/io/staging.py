"""Publish a downloaded cache file only once it is complete.

The ancillary caches are keyed by day (MERRA-2) or month (CAMEL), while the
callers work per scene -- an event's five scans of one day all ask for the
same file. Serially that is a cache hit; run the scans in parallel and it is
two processes downloading into one directory at the same instant, which used
to go wrong twice over:

* ``earthaccess.download`` names the granule after the remote object, so both
  processes wrote the same path, and the first to finish had its file deleted
  (``unlink``) or renamed out from under it by the second;
* the cache itself was written in place, so a third process could open a
  half-written netCDF -- silently, since the file exists and only its tail is
  missing.

``staged_download`` gives each call its own scratch directory inside the cache
directory and publishes the result with ``os.replace``. A cache path is then
either absent or complete, and duplicate downloads waste bandwidth instead of
corrupting each other.
"""

import os
import shutil
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def staged_download(cache: Path) -> Iterator[tuple[Path, Path]]:
    """Yield ``(scratch_dir, staged_path)`` for building ``cache``.

    Write the finished file to ``staged_path`` (and use ``scratch_dir`` for
    anything downloaded on the way); on a clean exit it is renamed onto
    ``cache``. The scratch directory is removed either way, so a failed or
    interrupted download leaves nothing behind but the missing cache. The
    rename is atomic because the scratch directory is a sibling of the cache,
    hence on its filesystem.
    """
    cache.parent.mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(dir=cache.parent, prefix=f".{cache.stem}-"))
    try:
        yield scratch, scratch / cache.name
        os.replace(scratch / cache.name, cache)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
