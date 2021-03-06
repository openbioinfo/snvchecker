.. SnvChecker documentation master file, created by
   sphinx-quickstart on Thu Oct 25 17:35:37 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SnvChecker's documentation!
======================================

Introduction
============

SnvChecker is an extension tool to call genotype for suscep sites, and will be firstly used in yigan pipeline & project.


Authors
=========
.. _authors:

    - denglh <denglonghui@venomics.com>
    - kongdeju <kongdeju@venomics.com>

Status
======

.. note::

    **not reviewed yet.**

Installation
============

use git to clone code::

    git clone git@192.168.1.251:/home/git/SnvChecker.git


Usage
=====

....
just type command::

    /path/to/locuschecker.py.py -h



must_args
---------

- method
    calling sambamba or freebayes

- bamfile
    path of sorted and indexed bam file

optinal args
------------

- variants
    path for risk information file
- reference
   path for the reference fasta
- prefix
   string for the sample id (name)
- outdir
   folder path for save output file (genetype output)


RUN
========

client way
-------

prepare requirement files and call script::

    /path/of/locuschecker.py freebayes -r <reference>   -v <variants>  [-o <output>]  [-p prefix] <bam>
    /path/of/locuschecker.py sambamba  [-r <reference>] -v <variants>  [-o <output>]  [-p prefix] <bam>



Tests
=========

check test report `here <http://192.168.1.4700:/dev-tests/SnvChecker/>`


Report
=========

check sample report `here <http://192.168.1.4700:/dev-report/SnvChecker/html/>`


Code
=======


.. toctree::
   :maxdepth: 1

   Guide <index>
   Code Docs <api/modules>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
