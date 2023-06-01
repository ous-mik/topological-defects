# Topological defects

## Manuscript

This repository contains Python scripts accompanying the following manuscript:

Emma Lång, Anna Lång, Pernille Blicher, Torbjørn Rognes, Paul Gunnar Dommersnes, and Stig Ove Bøe (2023)</br>
**Topology-guided polar ordering of epithelial monolayers** (submitted)

Please see the Methods section of the manuscript for details.


## Calculation of cell velocity

To calculate the velocity field using Particle Image Velocimetry (PIV), we employed a script that utilized the OpenPIV package. The average cell speed was extracted from the PIV data using the Python script `Defects.py`.

In addition, data on the average cell speed in monolayers subjected to stimulation, was extracted by single particle tracking carried out using the TrackMate plugin in ImageJ. The average cell speed per frame was calculated using the Python script `TrackMate_speed_data.py`.


## Detection and analysis of +1 and -1 topological defects

The +1 and -1 topological defects were detected using the Python script `Defects.py`. Dependencies: openpiv (version 0.24.2), numpy (version 1.23.5), matplotlib (version 3.7.1), pandas (version 1.5.3), scipy (version 1.10.1), imageio (version 2.26.1).


## Calculation of spatial correlation length

The correlation length was retrieved using the Python script `Spatial_correlation.py`.


## Calculation of vorticity

The Vorticity was retrieved using the Python script `Vorticity.py`.
