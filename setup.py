__author__ = 'dgaston'

"""Setup file and install script"""
import sys
import os
from setuptools import setup, find_packages

from ngsflow.version import __version__ as version


def write_version_py():
    version_py = os.path.join(os.path.dirname(__file__), 'bcbio', 'pipeline',
                              'version.py')
    try:
        import subprocess
        p = subprocess.Popen(["git", "rev-parse", "--short", "HEAD"],
                             stdout=subprocess.PIPE)
        githash = p.stdout.read().strip()
    except:
        githash = ""
    with open(version_py, "w") as out_handle:
        out_handle.write("\n".join(['__version__ = "%s"' % version,
                                    '__git_revision__ = "%s"' % githash]))

with open("requirements.txt", "r") as f:
    install_requires = [x.strip() for x in f.readlines() if not x.startswith(("bcbio-nextgen", "#"))]

# library-only install: enable skipping of scripts and requirements for conda builds
if "--record=/dev/null" in sys.argv:
    scripts = []
    install_requires = []
    zip_safe = True
else:
    zip_safe = False
    scripts = ['scripts/bcbio_nextgen.py', 'scripts/bcbio_setup_genome.py', 'scripts/bcbio_prepare_samples.py']

write_version_py()
setup(name="helenus",
      version=version,
      author="Dan Gaston",
      author_email="daniel.gaston@dal.ca",
      description="Pipeline infrastructure and services for processing next-generation sequencing data",
      long_description=(open('README.rst').read()),
      license="MIT",
      url="https://github.com/dgaston/nsha-mdx",
      packages=find_packages(),
      zip_safe=zip_safe,
      scripts=scripts,
      install_requires=install_requires)