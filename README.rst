========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |travis| |requires|
        |
    * - package
      - |version| |downloads| |wheel| |supported-versions| |supported-implementations|

.. |docs| image:: https://readthedocs.org/projects/ddb-ngsflow/badge/?style=flat
    :target: https://readthedocs.org/projects/ddb-ngsflow
    :alt: Documentation Status

.. |travis| image:: https://travis-ci.org/dgaston/ddb-ngsflow.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/dgaston/ddb-ngsflow

.. |requires| image:: https://requires.io/github/dgaston/ddb-ngsflow/requirements.svg?branch=master
    :alt: Requirements Status
    :target: https://requires.io/github/dgaston/ddb-ngsflow/requirements/?branch=master

.. |version| image:: https://img.shields.io/pypi/v/ddb-ngsflow.svg?style=flat
    :alt: PyPI Package latest release
    :target: https://pypi.python.org/pypi/ddb-ngsflow

.. |downloads| image:: https://img.shields.io/pypi/dm/ddb-ngsflow.svg?style=flat
    :alt: PyPI Package monthly downloads
    :target: https://pypi.python.org/pypi/ddb-ngsflow

.. |wheel| image:: https://img.shields.io/pypi/wheel/ddb-ngsflow.svg?style=flat
    :alt: PyPI Wheel
    :target: https://pypi.python.org/pypi/ddb-ngsflow

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/ddb-ngsflow.svg?style=flat
    :alt: Supported versions
    :target: https://pypi.python.org/pypi/ddb-ngsflow

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/ddb-ngsflow.svg?style=flat
    :alt: Supported implementations
    :target: https://pypi.python.org/pypi/ddb-ngsflow


.. end-badges

A toil based NGS workflow manager

* Free software: BSD license

Installation
============

::

    git clone https://github.com/dgaston/ddb-ngsflow
    cd ddb-ngsflow
    python setup.py install

Documentation
=============

https://ddb-ngsflow.readthedocs.org/

Development
===========

To run the all tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox
