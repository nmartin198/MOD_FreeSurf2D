# MOD_FreeSurf2D
**MOD_FreeSurf2D** is a computer program to solve the depth-averaged, shallow water equations in general situations. It was originally released in 2005 along with the article ["MOD_FreeSurf2D: a Matlab surface fluid flow code for rivers and streams"](https://doi.org/10.1016/j.cageo.2005.03.004).

In the original version, and the version provided here, **MOD_FreeSurf2D** is a set of [Matlab scripts](https://github.com/nmartin198/MOD_FreeSurf2D/tree/tsinout/src). These scripts require Matlab for execution. A release executable that is compiled, using the Matlab Compiler package, from the Matlab scripts is also available. It can be installed using the [Installer Package](https://github.com/nmartin198/MOD_FreeSurf2D/releases/tag/v2.0).

This branch differs from the main branch of the repository in that **MOD_FreeSurf2D** has been modified to use a [HDF5 file format](https://portal.hdfgroup.org/display/knowledge/What+is+HDF5). All input parameters are now contained within the input HDF5 file; all outputs are written to the same HDF5 file. The 

**MOD_FreeSurf2D** is designed to simulate water flow in rivers, streams, and shallow estuaries. Additional information on the required inputs, use of this program, and interpretation of results can be found in the [manual](https://github.com/nmartin198/MOD_FreeSurf2D/tree/tsinout/Docs).

## Requirements

For customizable and extendable use, **MOD_FreeSurf2D** requires Matlab. Running these scripts from within Matlab provides the user with unlimited access to customize, visualize, and analyze simulation results. The [Matlab scripts](https://github.com/nmartin198/MOD_FreeSurf2D/tree/tsinout/src) have been tested on Matlab 2020R.

A compiled executable is also available that can be installed using the [Installer Package](https://github.com/nmartin198/MOD_FreeSurf2D/releases/tag/v2.0). The Matlab Runtime is also installed by the Installer Package and allows execution of **MOD_FreeSurf2D** on computers, which do not have a Matlab license. The incorporation of the HDF5 file for inputs and outputs means that the user can now obtain results for the entire domain at specified time intervals and for every time step at specified point locations within the domain.

## Getting Started

First download the archive (and expand it) to a local hard drive. Add the MOD_FreeSurf2D directory to the current Matlab path. The program can be run from the Matlab command window by typing MOD_FreeSurf2D.

Three complete [examples](https://github.com/nmartin198/MOD_FreeSurf2D/tree/tsinout/Examples) are provided.

1. ConvS3
2. Reach1
3. DamBreak

Example [Jupyter Notebooks](https://jupyter.org/) are provided to illustrate creation of the HDF5 inputs with the formatting expected by **MOD_FreeSurf2D** and to provide an example of removing simulation results from an existing HDF5 file. Simulation results have to be removed from the HDF5 file before running a new simulation so that [Matlab](https://www.mathworks.com/products/matlab.html) does not throw an error when it attempts to overwrite existing information. These example Notebooks are available in the [JNotes](https://github.com/nmartin198/MOD_FreeSurf2D/tree/tsinout/Examples/JNotes) directory.

The total input and output for these examples is located in the [Examples directory](https://github.com/nmartin198/MOD_FreeSurf2D/tree/tsinout/Examples) under the appropriate subdirectory. In addition, plotting and output routines that have been customized for each scenario are also provided in the appropriate subdirectory. Each of these examples can be quickly run for testing purposes. To run **MOD_FreeSurf2D**, launch the executable from a command window. It will then prompt for the HDF5 file to use. Please see the main code branch for additional information on running the examples.

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

