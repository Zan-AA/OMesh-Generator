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

## Result
Here, we show two figures that resemble with the figures (Fig.1 and Fig. 2) shown in the paper.
![Fig1](https://github.com/Zan-AA/OMesh-Generator/blob/main/Fig1.png)

<p align="center">
  <img width="640" height="480" src="https://github.com/Zan-AA/OMesh-Generator/blob/main/Fig2.png">
</p>
