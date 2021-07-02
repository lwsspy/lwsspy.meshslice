MeshSlice
---------

Simple/small tool to plot and slice 3D meshes and create reproducible 2D slices.

+---------------+-------------------------------------------------------------------------------------+
| Documentation | https://lwsspy.github.io/lwsspy.meshslice/                                          |
+---------------+-------------------------------------------------------------------------------------+
| Deployment    | .. image:: https://img.shields.io/pypi/v/meshslice.svg?logo=python&logoColor=white) |
|               |     :target: https://pypi.org/project/meshslice/0.0.3/                              |
+---------------+-------------------------------------------------------------------------------------+
| Build Status  | .. image:: https://travis-ci.com/lwsspy/lwsspy.meshslice.svg?branch=main            |
|               |     :target: https://travis-ci.com/lwsspy/lwsspy.meshslice                          |
+---------------+-------------------------------------------------------------------------------------+
| License       | .. image:: https://img.shields.io/badge/License-GPLv3-blue.svg                      |
|               |     :target: https://www.gnu.org/licenses/gpl-3.0                                   |
+---------------+-------------------------------------------------------------------------------------+


Quick-Install
=============

.. code:: bash
    
    pip install lwsspy.meshslice



Quick-Usage
===========

.. code:: python
        
    from lwsspy.meshlice import MeshPlotSph, read_mesh
    M = read_mesh(<yourmeshfile>)
    Meshplot(M)






