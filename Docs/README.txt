README.txt file for MOD_FreeSurf2D
4 June 2004

1.  Introduction
2.  Requirements
3.  Installation
4.  Directory Structure
5.  Examples
6.  List of Script Files


1. Introduction

MOD_FreeSurf2D is a set of Matlab scripts that solve the
depth-averaged, shallow water equations in general situations.  These
scripts require Matlab for execution.  This program is designed to
simulate water flow in rivers, streams, and shallow estuaries.  The
program requires at least four input files to run: Depth.txt,
input.txt, Mann.txt, and Topo.txt.  Additional information on the
required input files and on the use of this program can be found in
the manual, MOD_FreeSurf2DManual.pdf, located in the Docs
directory. This code is described in "MOD_FreeSurf2D: a Matlab surface
fluid flow code for rivers and streams" by N. Martin and S. Gorelick.
Computers & Geosciences


2.  Requirements

MOD_FreeSurf2D should run on any working installation of Matlab 6.5
and higher, but we provide no warranty or guarantee of any kind.
Several of the Matlab functions used in the code are not included in
Matlab 5.x.  As a result, the scripts will not run on Matlab 5.x
without modification.


3.  Installation

a) Copy the MOD_FreeSurf2D directory (and subdirectories) to your hard
drive.  The files within the MOD_FreeSurf2D directory (not including
the files in the subdirectories) are the Matlab m-scripts that compose
the program.  As long as the MOD_FreeSurf2D directory is in the current
Matlab path, the program can be run from the Matlab command window by
typing MOD_FreeSurf2D.


4.  Directory Structure

MOD_FreeSurf2D directory is the main/top directory.

     Docs subdirectory holds documentation including:
                MOD_FreeSurf2DManual.pdf = program manual
                README.txt

     Examples subdirectory holds example simulation output for
          comparison to the output from simulations.  Includes a
          subdirectory for each example simulation:
                ConvS3/
                Reach1/
                DamBreak/

     sims subdirectory is to be used to run the three example
	  simulations.  Includes a subdirectory for each example
          simulation with the following files:
                ConvS3/
		     Depth.txt
                     Mann.txt
                     Topo.txt
                     TopoX.txt
                     TopoY.txt
                     cpLOT.m
                     fluid.m
                     input.txt
                Reach1/
		     Depth.txt
                     Mann.txt
		     Rch1_Comp_AVel.dat
		     Rch1_Comp_Depth.dat
                     Topo.txt
                     TopoX.txt
                     TopoY.txt
		     cOMPsTATgEN.m
                     cpLOT.m
                     fluid.m
                     input.txt
		     mODsTATgEN.m
                DamBreak/
		     Depth.txt
                     Mann.txt
                     Topo.txt
		     aFTERpLOT.m
		     dbinit.m
		     dbset.m
		     dbwrite.m
                     fluid.m
                     input.txt
		     loc1.txt
		     loc2.txt
		     loc3.txt
		     loc4.txt
		     loc5.txt
		     loc6.txt
		     loc8.txt
		     rEADiNPUT.m
		     sTATgEN.m


5.  Examples

Three complete examples are provided: ConvS3, Reach1, and DamBreak.
The total input and output for these examples is located in the
Examples directory under the appropriate subdirectory.  The input
files necessary to produce this output are located in the appropriate
subdirectory in the sims directory tree.  In addition, plotting and
output routines that have been customized for each scenario are also
provided in the appropriate sims subdirectory.  Consequently, each of
these examples can be quickly run (about 2 - 20 minutes) from the sims
subdirectory and the output checked against the output files in the
Examples subdirectory.  The specifics of running each example case are
listed below.

a) ConvS3 is a straight, sloping, rectangular channel case with flow
from left to right.  Water flows into the channel at the left (the
inflow boundary condition is fixed total water depth) and flows out of
the channel at the right (outflow boundary condition is radiation free
surface). Additional Matlab scripts are provided in this subdirectory
to provide a simple velocity contour plot at the end of the
simulation.  To run this example, 1) start Matlab 2) set the current
directory to be the sims/ConvS3 directory 3) add the MOD_FreeSurf2D
directory to the current Matlab path 4) copy the file fluid.m from the
ConvS3 subdirectory to the MOD_FreeSurf2D directory, and 5) type
MOD_FreeSurf2D at the Matlab command line.  Copying the file fluid.m
will overwrite the version of fluid.m that originally exists in the
MOD_FreeSurf2D directory.  A copy of the original file is stored in
the MOD_FreeSurf2D directory as org_fluid.m.

b) Reach1 comes from a USGS study of part of the Kootenai River, ID.
This is an approximately 500 meter reach.  Again, flow is from left to
right.  Inflow boundary conditions are specified flux.  Outflow
boundary conditions are radiation free surface.  The simulation in the
Reach1 directory employs a spatial discretization of 10 meters and a
time step of 15.0 seconds.  Additional Matlab scripts are included in
the Reach1 directory to provide plots of the simulation and an error
analysis of the simulated results.  The results are compared to the
data in the files Rch1_Comp_AVel.dat and Rch1_Comp_Depth.dat for the
error analysis.  To run this example, follow steps 1-5 in the ConvS3
example except that need to set the current directory to be
sims/Reach1.

c) DamBreak comes from the dam-break style flume experiment of Bellos
et al. (1992).  The flume in this experiment is about 21 meters long
and about 1.4 meters wide.  The flume has a slope of 0.002.  At the
start of the experiment, the flume below the dam is dry.  This example
simulation lasts for about 70.0 seconds of "real time".  The files
loc1.txt,...,loc8.txt contain the measured depths at the seven
locations for which Bellos et al. (1992) provide water depth data.
Again, extra scripts are included to do some plotting and simulation
accuracy analysis.  To run this example follow steps 1-3 in the ConvS3
example (but set the current directory to be sims/DamBreak).  As in
the other examples the fourth step is to copy the fluid.m file in the
sims/DamBreak directory to the MOD_FreeSurf2D directory.  The
rEADiNPUT.m file also needs to be copied from the sims/DamBreak
directory to the MOD_FreeSurf2D directory.  A copy of the original
rEADiNPUT.m file is stored in the MOD_FreeSurf2D directory as
org_rEADiNPUT.m.  Now to run this example, type MOD_FreeSurf2D at the
Matlab command prompt.


6.  List of script files

MOD_FreeSurf2D.m - main wrapper that calls all other scripts.

aLLOCaDJvOL.m - creates vectors holding variable values for adjacent
	  computational volumes.
aLLOCfACE.m - creates vectors holding variables values corresponding
	  to a particular face of each computational volume.
aVEdEPTHcALC.m - calculates the average water depth in each
	  computational volume.
cHEZYcALC.m - calculates the Chezy friction coefficients.
cOLaDJ.m - ensures the column index lies inside of domain.
cOLbaDJ.m - ensures the column index lies inside of domain.
dELTAcALC.m - calculates part of the RHS of the free surface system of
	  equations.
dOMsYSbOUND.m - enforces boundary conditions for wet and dry domain
	  boundaries.
eULiNT.m - does path line tracing with Euler integration.
fLUIDiNIT.m - initializes global variables.
fLUIDmASSsETUP.m - sets up Mass.txt output file.
fLUIDmASSwRITE.m - writes values to Mass.txt.
fREEsURF.m - wrapper for free surface calculations.
fUXcALC.m - wrapper for semi-Lagrangian advection operator
	  calculations at x-faces.
fVYcALC.m - wrapper for semi-Lagrangian advection operator
	  calculations at y-faces.
fcORcALC.m - calculates Coriolis f-term with f-plane model.
fluid.m - "main" fluid flow program.
fsbOUNDaDJ.m - adjusts free surface values for dry and wet domain
	  boundary locations.
futERMcALC.m - calculates Coriolis term for x-face semi-Lagrangian
	  advection operator.
fvtERMcALC.m - calculates Coriolis term for y-face semi-Lagrangian
	  advection operator.
gENfACEiNDEX.m - generates indexes for for face and adjacent volume
	  vectors.
iNTbOUND.m - adjusts calculated particle locations (semi-Lagrangian
	  advection representation) for boundaries.
iNdOMAINmASS.m - calculates mass of water in simulation domain.
lOCATIONcALC.m - calculates particle volume index location.
lhscALC.m - generates the LHS of the free surface system of equations.
mASSfLUXcALC.m - determines the net mass of water entering the
	  simulation domain.
mainfluid.m - wrapper program that calls the fluid flow calculations
	  completed in every time step.
nORMd.m - normalizes free surface system of equations.
org_fluid.m - 'original' version of fluid.m
org_rEADiNPUT.m - 'original' version of rEADiNPUT.m.
qiNfLUXbc.m - handles Dirichlet influx boundary conditions.
rADfLUXbc.m - handles radiation flux boundary conditions.
rADoRLfSbc.m - handles radiation free surface boundary conditions.
rADoRLsET.m - set up program to deal with radiation free surface
	  boundary conditions.
rADvELbc.m - handles radiation velocity boundary conditions.
rADvELsET.m - set up program to deal with radiation velocity boundary
	  conditions.
rEADiNPUT.m - reads the input file, input.txt.
rIToUTPUT.m - writes selected output to text files.
rOWaDJ.m - ensures that row indexes are inside simulation domain.
rOWbaDJ.m - ensures that row indexes are inside simulation domain.
rhscALC.m - generates the RHS of the free surface system of equations.
rkiNT.m - tracks particle path lines with Runge-Kutta integration.
sEMIapATH.m - tracks particle path lines with the semi-analytic method
	  of Pollock (1988).
sETafACES.m - generates 'A' calculation vectors for volume faces.
sETdOMAIN.m - creates the domain layout for a simulation.
sETgfACES.m - generates 'G' calculation vectors for volume faces.
sETiNITIALvECTORS.m - sets some initial vector quantities from the
          domain layout.
sETpRECISION.m - sets the smallest numbers that can be employed in
          calculations within the program.
sETtOPO.m - calculates volume face topographic elevations.
sETtOTALdEPTH.m - updates total water depth.
tOTmANNINGcALC.m - prepares Manning's n values for calculations.
tdEPdIRCbc.m - enforces Dirichlet total depth boundary conditions.
tdEPsYSbc.m - sets Dirichlet total depth boundary conditions in free
          surface system of equations.
ucUBIClATiNT.m - does cubic Lagrange polynomial interpolation for
          x-face values.
udIFFcALC.m - calculates eddy diffusivity contribution to x-face
          advective operator.
vELcALC.m - updates velocity values.
vELdIRCbc.m - enforces Dirichlet velocity domain boundaries.
vELpARAMsET.m - sets velocity related calculation vectors.
vELsET.m - sets initial velocity if specified with BU.txt and BV.txt.
vOLUMEfIND.m - finds the current computational volume index for each
          tracked particle.
vOLtOXfACE.m - upwinds volume center values to volume x-faces.
vOLtOYfACE.m - upwinds volume center values to volume y-faces.
vOLtOfACE.m - wrapper that calls vOLtOXfACE.m and vOLtOYfACE.m
vcUBIClATiNT.m - does cubic Lagrange polynomial interpolation for
          y-face values.
vdIFFcALC.m - calculates eddy diffusivity contribution to y-face
          advective operator.
xINDEXfIND.m - finds x-face indexes for bilinear interpolation.
yINDEXfIND.m - finds y-face indexes for bilinear interpolation.
