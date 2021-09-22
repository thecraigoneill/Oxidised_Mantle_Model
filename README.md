# Oxidised_Mantle_Model
Rheological model for Aspect - combining homologous scaling for bridgmanite and magnesiowuestite. From O'Neill and Aulbach, 2021.

<a href="https://zenodo.org/badge/latestdoi/409012910"><img src="https://zenodo.org/badge/409012910.svg" alt="DOI"></a>

![Brig14_2025bMyr](https://user-images.githubusercontent.com/30849698/134265601-7d9b3eb3-5c03-42eb-a313-8343317c4426.png)

This module is designed to work with Community mantle convection solver Aspect:

https://github.com/geodynamics/aspect

For the Aspect manual, and instructions on compiling a module, see the manual: www.math.clemson.edu/~heister/manual.pdf
(see pg. 280 on for details on plugins and compilation).

Once you have a working Aspect (+ deal.II etc) distribution, go into the subfolder with the cmake/cc/h files, and type:

cmake -DAspect_DIR=/path/to/aspect/build/ .

make

This should make the source file (so), and link to an aspect executable. The supplied prm input file has an example implementation, and shows how to call the plugin.

Also supplied is an example pbs for submission to cluster environments - note this is machine dependent, obviously. 
