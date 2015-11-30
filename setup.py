__author__ = 'dgaston'

"""Setup file and install script"""
import os
import io
from setuptools import setup, find_packages

from ngsflow.version import __version__ as version

here = os.path.abspath(os.path.dirname(__file__))


def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)


long_description = read('README.md')

with open("requirements.txt", "r") as f:
    install_requires = [x.strip() for x in f.readlines() if not x.startswith(("ddbio-", "#"))]

setup(name="ddbio-ngsflow",
      version=version,
      author="Dan Gaston",
      author_email="daniel.gaston@gmail.com",
      description="Pipeline infrastructure and services for processing next-generation sequencing data",
      long_description=long_description,
      license="MIT",
      url="https://github.com/dgaston/ddbio-ngsflow",
      packages=find_packages(),
      install_requires=install_requires)
