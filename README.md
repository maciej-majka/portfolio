# Simulation engine for polymeric chain driven by Spatially Correlated Noise (SCN)
## Description
This is a simulation engine for 2D polymeric chain, driven by SCN. Initially created for my MSc thesis and further developed for my early PhD, I worked on it somwhere bewteen 2010 and 2012. It was used in my research on the role of SCN in physical systems. The results generated with this code were discussed in a few of my publications[^1][^2][^3][^4], where you can find details about the physical and mathematical context.

_Warning: This old code is uploaded 'as is'. It is poorely commented (mostly in Polish), full of abandoned functions and commented-out sections, probably an anti-example of how serviceable code should be written... But it worked well and served its purpose in its time :)_

The code simulates the following physical situation:
<p align="center">
<img src="polymer.png" width="200">
</p>
Here, the red line is our polymeric chain and the arrows in the background are the snapshot of SCN vectors. Notice, that vectors close to each other (within certain 'correlation length') are 'similar' in direction and length. This is different from the non-correlated noise, for which the relative direction of two nearby vectors would be completely random. In order to simulate this system, the following components had to be implemented into to the code:

* Structures describing the polymer state and physical variables
* Methods for integrating the equations of motion
* Matrix operations for generating SCN
* Animations
* Interface for outputing data

## Structures describing the polymer state and physical variables

File set.h contains the declarations of objects consituting the simulated system and functions realizing the simulation.

* structure bead: It contains the state of a single polymer's bead e.g. its position, interaction constants etc.
* structure polymer: It describes the whole polymeric chain. Its most important member is a pointer under which a list of beads is allocated. It also contains data such as position of chain mass center etc.
* class set: This is the representation of the entire system and a set of functions actually performing the simulations. Important members include:
  - Pointer '*polymers' under which a list of polymers is allocated.
  - Function 'solution(...)', which performs integration of equations of motion over given time period and outputs the data
  - Pointers '*matrix' and '*matrixLT' where dynamic correlation matrix and its Cholesky decompositions are allocated   

## Methods for integrating the equation of motion

The simulation of physical system requires the numerical solution of differential equations. Since Stochastic Differential Equations must be treatet with special algorithms, it was necessary to implement this procedure 'by hand'. Initially, the stochastic Verlet method was applied (whose remanants still exist in class set), but eventually the stoachstic Runge-Kutta 4th order method was implemented as more reliable for the relatively complex interactions contained in the model. This method is called within function set::solution.

## Matrix operations for generating SCN
In order to generate the vector of SCN one has to multiply a vector of non-correlated noise by the matrix resulting from Cholesky decomposition of correlation matrix. Class matrix and its utilities were written mostly for training purposes. Non-correlated noise was generated using Mersenne-Twister generator from GSL library.

## Animations
Animations allowed visual instepction of the code performance in single test runs. Two variants were written:

* plot_anim.h and anim_fun.cpp : This version required Gnuplot. It transferred current polymer position to external file and called Gnuplot to display it. Then updated the positions and send refresh prompt to Gnuplot
* animGL.h and animGL.cpp : This is OpenGL-based version (using freeglut library). The main simulation programme called separate, auxilary program for displaying animation and periodically outputted current polymer positions to external file. The auxilarly program, containing OpenGL loop, read the external file to update the displayed positions of the polymer beads.

The function feeding data to external file can be found at the end of set_tools.cpp.

# File description

* animGL.h, animGL.cpp: function for OpenGL-based animations
* anim_fun.cpp, plot_anim.h: functions for Gnuplot-based animtions
* forces.cpp: functions calculating energies and forces acting on polymer chain
* functions.cpp: defines noise correaltion function and various utility functions
* main.cpp: main file
* makefile: makefile used for compilations and removal of auxilarly files
* matrix.h, matrix_LT.cpp, matrix_ctors.cpp, matrix_fun.cpp: definition and implementation of matrix class. matrix_LT contains code for Cholesky decomposition, used for SCN generation
* set.h, set_ctors.cpp: define class set
* set_statistic.cpp: defines functions calculating statistical data from simulations
* set_tools.cpp: defines functions for outputing simulation data

[^1]: [M. Majka, P. F. Góra, ‘Polymer shape dynamics induced by spatially correlated noise’, Acta Phys. Pol. B, 43, 5, 1133 (2012)](http://web.a.ebscohost.com/abstract?direct=true&profile=ehost&scope=site&authtype=crawler&jrnl=05874254&AN=76482102&h=vQ3WB8hyb8pdYjbzvsWv6OdyFnl8cFEkHO4%2fiukVCmHp2kiWxuROMt1ic8rpbZa3Q6BHTSnrsCmEe0WrQGIN%2bw%3d%3d&crl=c&resultNs=AdminWebAuth&resultLocal=ErrCrlNotAuth&crlhashurl=login.aspx%3fdirect%3dtrue%26profile%3dehost%26scope%3dsite%26authtype%3dcrawler%26jrnl%3d05874254%26AN%3d76482102)
[^2]: [M. Majka, P. F. Góra, ‘Polymer unfolding and motion synchronization induced by spatially correlated noise’, Phys. Rev. E,  86, 5, 051122 (2012)](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.86.051122)
[^3]: [M. Majka, P. F. Góra, ‘Reinterpreting polymer unfolding effect induced by spatially correlated noise’, Acta Phys. Pol. B, 44, 5, 1099 (2013)](http://web.b.ebscohost.com/abstract?direct=true&profile=ehost&scope=site&authtype=crawler&jrnl=05874254&AN=88950514&h=7nct0WEPBJizKrO%2bYcidZI9vaBcNfhMMJPYjkNehyxDOkct7sWinj24GCBrPEpSJoGsvgW5%2bvsfLYOr4WYFhxA%3d%3d&crl=c&resultNs=AdminWebAuth&resultLocal=ErrCrlNotAuth&crlhashurl=login.aspx%3fdirect%3dtrue%26profile%3dehost%26scope%3dsite%26authtype%3dcrawler%26jrnl%3d05874254%26AN%3d88950514)
[^4]: [M. Majka, P. F. Góra, ‘Non-Gaussian polymers described by alpha-stable chain statistics: model, effective interactions in binary mixtures and application to on-surface separation’, Phys. Rev. E, 91, 5, 052602 (2015)](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.052602)
