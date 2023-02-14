# OMesh Generator
## Desciption
This script details the generation of 2D Euler-type structured O-meshes using Karman–Trefftz conformal transformation. Details of the generation procedure as well as the theoretical background can be found in the paper by Vassberg & Jameson:

```
@article{vassberg2010pursuit,
  title={In pursuit of grid convergence for two-dimensional Euler solutions},
  author={Vassberg, John C and Jameson, Antony},
  journal={Journal of Aircraft},
  volume={47},
  number={4},
  pages={1152--1166},
  year={2010}
}
```

## Installation
This script requires the minimum python3 and matplotlib functionalities. For saving the mesh, the script provides one option of using plot3d (https://github.com/nasa/plot3d_utilities) which can be installed via
```
pip3 install plot3d
```
The user can also modify the script to generate other meshes (such as using meshio).

## Example
The simplest command is to use the default setting
```
python3 Omesh.py 
```
which generates a 32 x 32 O-mesh. The final result is also shown in the form of a plot. 

The script takes in command line arguments such as `NC` (Number of cells in both radial and circumferential directions). The user can also mute the plot by using `-p 0`. E.g., the following command generates a 64 x 64 O-mesh, without showing the final plot
```
python3 Omesh.py -NC 64 -p 0
```
To save the final result into a plot3d mesh file (require the plot3d python module), `-s 1` command is required
```
python3 Omesh.py -NC 64 -s 1
```
Users are directed to the `CommandLineParser` class in the script for details of the input arguments.

## Result
Here, we show two figures that resemble with the figures (Fig.1 and Fig. 2) shown in the paper.
![Fig1](https://github.com/Zan-AA/OMesh-Generator/blob/main/Fig1.png)

### A 32x32 O-Mesh with close-up view near the airfoil.
<p align="center">
  <img width="640" height="480" src="https://github.com/Zan-AA/OMesh-Generator/blob/main/Fig2.png">
</p>
