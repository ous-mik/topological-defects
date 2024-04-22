# Topological defects

## Publication

This repository contains Python scripts accompanying the following paper:

Emma Lång, Anna Lång, Pernille Blicher, Torbjørn Rognes, Paul Gunnar Dommersnes, and Stig Ove Bøe (2024)</br>
**Topology-guided polar ordering of collective cell migration**</br>
Science Advances, 10 (16), eadk4825</br>
doi: [10.1126/sciadv.adk4825](https://doi.org/10.1126/sciadv.adk4825)

Please see the Methods section of the paper for details.


## Scripts

All Python scripts are located in the folder `scripts`. The paths to the input and output files will need to be modified within the scripts.


### Calculation of cell velocity

To calculate the velocity field using Particle Image Velocimetry (PIV), we employed a script that utilized the OpenPIV package. The average cell speed was extracted from the PIV data using the Python script `Defects.py`.

In addition, data on the average cell speed in monolayers subjected to stimulation, was extracted by single particle tracking carried out using the TrackMate plugin in ImageJ. The average cell speed per frame was calculated using the Python script `TrackMate_speed_data.py`.


### Detection and analysis of +1 and -1 topological defects

The +1 and -1 topological defects were detected using the Python script `Defects.py`. Dependencies: openpiv (version 0.24.2), numpy (version 1.23.5), matplotlib (version 3.7.1), pandas (version 1.5.3), scipy (version 1.10.1), and imageio (version 2.26.1).


### Calculation of spatial correlation length

The correlation length was retrieved using the Python script `Spatial_correlation.py`.


### Calculation of vorticity

The Vorticity was retrieved using the Python script `Vorticity.py`.


### Plotting a normalized tensor field and a velocity field overlaid on a microscopy image

These plots were created using the Python script `tensor_overlay3.py`.


### Plotting the angles between the velocity and orientation of cells

These plots were created using the Python script `Tensor_velocity_angles2_batch.py`.
