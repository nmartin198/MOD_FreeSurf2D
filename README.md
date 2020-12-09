# MOD_FreeSurf2D
**MOD_FreeSurf2D** is a computer program to solve the depth-averaged, shallow water equations in general situations. It was originally released in 2005 along with the article ["MOD_FreeSurf2D: a Matlab surface fluid flow code for rivers and streams"](https://doi.org/10.1016/j.cageo.2005.03.004).

In the original version, and the version provided here, **MOD_FreeSurf2D** is a set of [Matlab scripts](https://github.com/nmartin198/MOD_FreeSurf2D/tree/main/src). These scripts require Matlab for execution. A release executable that is compiled, using the Matlab Compiler package, from the Matlab scripts is also available. It can be installed using the [Installer Package](https://github.com/nmartin198/MOD_FreeSurf2D/tree/main/Release).

**MOD_FreeSurf2D** is designed to simulate water flow in rivers, streams, and shallow estuaries. It requires at least four input files to run.

1. Depth.txt (initial water depth across the domain)
2. input.txt (input parameter specification)
3. Mann.txt (Manning's n for each computational grid cell)
4. Topo.txt (ground surface elevation for each grid cell)

Additional information on the required input files and on the use of this program can be found in the [manual](https://github.com/nmartin198/MOD_FreeSurf2D/tree/main/Docs).

## Requirements

For customizable and extendable use, **MOD_FreeSurf2D** requires Matlab. Running these scripts from within Matlab provides the user with unlimited access to customize, visualize, and analyze simulation results. The [Matlab scripts](https://github.com/nmartin198/MOD_FreeSurf2D/tree/main/src) have been tested on Matlab 2020R.

A compiled executable is also available that can be installed using the [Installer Package](https://github.com/nmartin198/MOD_FreeSurf2D/tree/main/Release). The Matlab Runtime is also installed by the Installer Package and allows execution of **MOD_FreeSurf2D** on computers, which do not have a Matlab license. The drawback to the compiled executable is that the user cannot customize the program outputs using the Matlab visualization and analysis tools.

## Getting Started

First download the archive (and expand it) to a local hard drive. Add the MOD_FreeSurf2D directory to the current Matlab path. The program can be run from the Matlab command window by typing MOD_FreeSurf2D.

Three complete [examples](https://github.com/nmartin198/MOD_FreeSurf2D/tree/main/Examples) are provided.

1. ConvS3
2. Reach1
3. DamBreak

The total input and output for these examples is located in the [Examples directory](https://github.com/nmartin198/MOD_FreeSurf2D/tree/main/Examples) under the appropriate subdirectory. In addition, plotting and output routines that have been customized for each scenario are also provided in the appropriate sims subdirectory. Each of these examples can be quickly run for testing purposes. The specifics of running each example case are listed below.

### ConvS3

ConvS3 is a straight, sloping, rectangular channel case with flow from left to right.  Water flows into the channel at the left (the inflow boundary condition is fixed total water depth) and flows out of the channel at the right (outflow boundary condition is radiation free surface). Additional Matlab scripts are provided in this subdirectory
to generate a simple velocity contour plot at the end of the simulation. To run this example:

1. Start Matlab
2. Set the current directory to be the Examples/ConvS3 directory
3. Type MOD_FreeSurf2D at the Matlab command line.

### Reach1

Reach1 comes from a USGS study of part of the Kootenai River, ID. This is an approximately 500 meter reach.  Again, flow is from left to right.  Inflow boundary conditions are specified flux.  Outflow boundary conditions are radiation free surface.  The simulation in the Reach1 directory employs a spatial discretization of 10 meters and a time step of 15.0 seconds.  Additional Matlab scripts are included in the Reach1 directory to provide plots of the simulation and an error analysis of the simulated results. Simulated results are compared to the data in the files Rch1_Comp_AVel.dat and Rch1_Comp_Depth.dat for the error analysis.

To run this example:

1. Start Matlab
2. Set the current directory to be the Examples/Reach1 directory
3. Type MOD_FreeSurf2D at the Matlab command line.

### Dambreak

DamBreak comes from a dam-break style flume experiment. The flume in this experiment is about 21 meters long and about 1.4 meters wide. The flume has a slope of 0.002. At the start of the experiment, the flume below the dam is dry. This example simulation lasts for about 70.0 seconds of "real time". The files loc1.txt,...,loc8.txt contain the measured depths at seven locations for which experimental water depth data are available. Again, extra scripts are included to do some plotting and simulation accuracy analysis.

To run this example: 

1. Start Matlab
2. Set the current directory to be the Examples/Dambreak directory
3. Type MOD_FreeSurf2D at the Matlab command line.


## Contributing

The author is happy to accept contributions to the project. It might be easiest to contact me prior to starting your contribution in case there are any initial suggestions or direction that I can provide.

The general procedure for contributing is as follows.

- Fork the project
- Make your changes
- Submit a pull request
    - It is important to have a conversation when opening a pull request. Describe your change and why it should be accepted.

## Author

* **Nick Martin** nick.martin@stanfordalumni.org

## License

This project is licensed under the GNU Affero General Public License v.3.0 - see the [LICENSE](LICENSE) file for details.

