function dbwrite
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% dbwrite writes the water depth output to files and plots these values
% over time.   This script is for simulation of the dambreak style flume
% experiments of Bellos et al. (1992).

global DepthLoc1 DepthLoc2 DepthLoc3 DepthLoc4 DepthLoc5 DepthLoc6 DepthLoc7
global DepthLoc8 TimeCounter

f1 = 'DL1.txt';
f2 = 'DL2.txt';
f3 = 'DL3.txt';
f4 = 'DL4.txt';
f5 = 'DL5.txt';
f6 = 'DL6.txt';
f7 = 'DL7.txt';
f8 = 'DL8.txt';
f9 = 'TC.txt';
save(f1,'DepthLoc1','-ASCII');
save(f2,'DepthLoc2','-ASCII');
save(f3,'DepthLoc3','-ASCII');
save(f4,'DepthLoc4','-ASCII');
save(f5,'DepthLoc5','-ASCII');
save(f6,'DepthLoc6','-ASCII');
save(f7,'DepthLoc7','-ASCII');
save(f8,'DepthLoc8','-ASCII');
save(f9,'TimeCounter','-ASCII');

return;
%EOF