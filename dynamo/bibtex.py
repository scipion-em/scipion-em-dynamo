# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     you (you@yourinstitution.email)
# *
# * your institution
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
@article{CASTANODIEZ2012139,
title = "Dynamo: A flexible, user-friendly development tool for subtomogram averaging of cryo-EM data in high-performance computing environments",
journal = "Journal of Structural Biology",
volume = "178",
number = "2",
pages = "139 - 151",
year = "2012",
note = "Special Issue: Electron Tomography",
issn = "1047-8477",
doi = "https://dx.doi.org/10.1016/j.jsb.2011.12.017",
url = "http://www.sciencedirect.com/science/article/pii/S1047847711003650",
author = "Daniel Castaño-Díez and Mikhail Kudryashev and Marcel Arheit and Henning Stahlberg",
keywords = "Subtomogram averaging, Single Particle Tomography, High-performance computing, GPU computing, Classification",
abstract = "Dynamo is a new software package for subtomogram averaging of cryo Electron Tomography (cryo-ET) data with three main goals: first, Dynamo allows user-transparent adaptation to a variety of high-performance computing platforms such as GPUs or CPU clusters. Second, Dynamo implements user-friendliness through GUI interfaces and scripting resources. Third, Dynamo offers user-flexibility through a plugin API. Besides the alignment and averaging procedures, Dynamo includes native tools for visualization and analysis of results and data, as well as support for third party visualization software, such as Chimera UCSF or EMAN2. As a demonstration of these functionalities, we studied bacterial flagellar motors and showed automatically detected classes with absent and present C-rings. Subtomogram averaging is a common task in current cryo-ET pipelines, which requires extensive computational resources and follows a well-established workflow. However, due to the data diversity, many existing packages offer slight variations of the same algorithm to improve results. One of the main purposes behind Dynamo is to provide explicit tools to allow the user the insertion of custom designed procedures – or plugins – to replace or complement the native algorithms in the different steps of the processing pipeline for subtomogram averaging without the burden of handling parallelization. Custom scripts that implement new approaches devised by the user are integrated into the Dynamo data management system, so that they can be controlled by the GUI or the scripting capacities. Dynamo executables do not require licenses for third party commercial software. Sources, executables and documentation are freely distributed on http://www.dynamo-em.org."
}

"""
