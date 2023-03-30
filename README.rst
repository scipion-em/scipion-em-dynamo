=========================
Scipion plugin for Dynamo
=========================

.. image:: https://img.shields.io/pypi/v/scipion-em-dynamo.svg
        :target: https://pypi.python.org/pypi/scipion-em-dynamo
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-dynamo.svg
        :target: https://pypi.python.org/pypi/scipion-em-dynamo
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-dynamo.svg
        :target: https://pypi.python.org/pypi/scipion-em-dynamo
        :alt: Supported Python versions

.. image:: https://img.shields.io/pypi/dm/scipion-em-dynamo
        :target: https://pypi.python.org/pypi/scipion-em-dynamo
        :alt: Downloads

This plugin allows to use Dynamo_ - an environment for subtomogram averaging cryo-ET data - within the Scipion framework.

============
Installation
============
The plugin can be installed in user (stable) or developer (latest, may be unstable) mode:

**1. User (stable) version:**:

.. code-block::

    scipion3 installp -p scipion-em-dynamo

**2. Developer (latest, may be unstable) version:**:

* Clone the source code repository:

.. code-block::

    git clone https://github.com/scipion-em/scipion-em-dynamo.git

* Move to devel branch:

.. code-block::

    git checkout devel

* Install:

.. code-block::

    scipion3 installp -p local/path/to/scipion-em-dynamo --devel

=========
Protocols
=========
The integrated protocols are:

1. bin tomograms: Reduce the size of a set of tomograms by a binning factor.

2. vectorial picking: Manual vectorial picking using Dynamo dtmslice GUI.

3. subtomogram extraction: Extraction of subtomograms.

4. import subtomograms from Dynamo .tbl files.

5. model workflow: apply a geometry model to a set of meshes, normally the result of the vectorial picking protocol.

6. subtomogram alignment.

=====
Tests
=====

The installation can be checked out running some tests. To list all of them, execute:

.. code-block::

     scipion3 tests --grep dynamo

To run all of them, execute:

.. code-block::

     scipion3 tests --grep dynamo --run

To run a specific test, for example, the tests to check the protocol for averaging subtomograms (the following command
can be copied from the test list displayed when listing the tests, as explained above):

.. code-block::

    scipion3 tests dynamo.tests.test_dynamo_average_subtomos.TestDynamoAverageSubtomograms

========
Tutorial
========

A tutorial with examples about how to use Dynamo within Scipion can be found here_.

==========
References
==========

* `Dynamo: A flexible, user-friendly development tool for subtomogram averaging of cryo-EM data in high-performance computing environments. <https://dx.doi.org/10.1016/j.jsb.2011.12.017>`_
  Daniel Castaño-Díez et al., Journal of Structural Biology, 2012.


===================
Contact information
===================

If you experiment any problem, please contact us here: scipion-users@lists.sourceforge.net or open an issue_.

We'll be pleased to help.

*Scipion Team*

.. _Dynamo: https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Main_Page
.. _here: https://scipion-em.github.io/docs/release-3.0.0/docs/user/tutorials/tomo/Picking_tutorial_lite/dynamo-tutorial-picking-lite.html#tutorial-picking
.. _issue: https://github.com/scipion-em/scipion-em-dynamo/issues
