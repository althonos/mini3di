#!/usr/bin/env python
# coding: utf-8

import configparser
import setuptools
from setuptools.command.sdist import sdist as _sdist

try:
    import astor
except ImportError as err:
    astor = err


class sdist(_sdist):
    """A `sdist` that generates a `pyproject.toml` on the fly.
    """

    def run(self):
        # build `pyproject.toml` from `setup.cfg`
        c = configparser.ConfigParser()
        c.add_section("build-system")
        c.set("build-system", "requires", str(self.distribution.setup_requires))
        c.set("build-system", 'build-backend', '"setuptools.build_meta"')
        with open("pyproject.toml", "w") as pyproject:
            c.write(pyproject)
        # run the rest of the packaging
        _sdist.run(self)


if __name__ == "__main__":
    setuptools.setup(cmdclass=dict(sdist=sdist))
