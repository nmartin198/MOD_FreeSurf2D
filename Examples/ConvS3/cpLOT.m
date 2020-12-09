function cpLOT
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% cpLOT provides a planview contour plot of calculated velocities for
% a stream reach.

global DIRD DX DY ENDTIME FILENAME FLUID_DT H INFO MN NUMCOLS NUMNODES
global NUMROWS PREC STARTTIME THETA TOPONAME TIMESTEP UAve VAve
global XINDEX YINDEX

% local variables.
NVel = zeros(NUMNODES,1);
CalVel = zeros(NUMROWS,NUMCOLS);
MVel = zeros(NUMROWS,NUMCOLS);

xloc = XINDEX(1)+0.5*DX:DX:XINDEX(NUMCOLS)+0.5*DX;
yloc = YINDEX(1)+0.5*DX:DX:YINDEX(NUMROWS)+0.5*DY;
[X,Y] = meshgrid(xloc,yloc);

NVel = sqrt(UAve.^2 + VAve.^2);
for i=1:NUMROWS
   for j=1:NUMCOLS
      CalVel(i,j) = NVel(((i-1)*NUMCOLS)+j);
   end
end

v = 0.1:0.1:2.0;

figure(2);
hold on;
[D,i] = contour(X,Y,CalVel,v);
% NM modification 2020-12-8
e = clabel(D);
%e = clabel(D,i);
caxis([v(1) v(20)]);
title('ConvS3 Example Simulation');
xlabel('Easting in meters','FontSize',12);
ylabel('Northing in meters','FontSize',12);
axis([XINDEX(1) XINDEX(NUMCOLS+1) YINDEX(1) YINDEX(NUMROWS+1)]);
%axis square;
%set(gca,'XTick',0:10:250);
%set(gca,'YTick',0:10:250);
%set(gca,'XTickLabel',int2str((0:10:250)'),'FontSize',8);
%set(gca,'YTickLabel',int2str((0:10:250)'),'FontSize',8);
%set(gca,'XGrid','on','GridLineStyle','-');
%set(gca,'YGrid','on','GridLineStyle','-');
grid on;
saveas(gcf,'ConvS3Example.fig','fig');
hold off;

clear D CalVel e NVel X xloc Y yloc;
return;
%EOF
