function cpLOT
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% cpLOT provides a planview contour plot of calculated velocities for
% a stream reach.

global DX DY H INFO MN NUMCOLS NUMNODES
global NUMROWS NUMINCX NUMINCY PREC u v THETA
global XINDEX YINDEX ROWBEGIN ROWEND XINC

% local variables.
u1 = double(zeros(NUMNODES,1));
u3 = double(zeros(NUMNODES,1));
v2 = double(zeros(NUMNODES,1));
v4 = double(zeros(NUMNODES,1));
UAv = double(zeros(NUMNODES,1));
VAv = double(zeros(NUMNODES,1));
NVel = double(zeros(NUMNODES,1));
CalVel = double(zeros(NUMROWS,NUMCOLS));
CalDep = double(zeros(NUMROWS,NUMCOLS));
MVel = double(zeros(NUMROWS,NUMCOLS));
PUv = double(zeros(NUMROWS,NUMCOLS+1));
PVv = double(zeros(NUMROWS+1,NUMCOLS));

fnv = 'AVelxyz.dat'; % average velocity i,j,k file
fnd = 'Depxyz.dat'; % average depth i,j,k file
fpeps = 'Plots.eps';
fpmat = 'Plots.fig';
f1 = fopen(fnv,'w+');
f2 = fopen(fnd,'w+');

XMin = 543380;
XMax = 543890;
YMin = 5403660;
YMax = 5403990;
PlotX = 543400:100:543800;
PlotY = 5403700:50:5403950;

xloc = (XMin:DX:(XMin + XINDEX(NUMCOLS)));
x_xloc = ((XMin - 0.5*DX):DX:((XMin + 0.5*DX) + XINDEX(NUMCOLS)));
yloc = (YMin:DY:(YMin + YINDEX(NUMROWS)));
y_yloc = ((YMin - 0.5*DY):DY:((YMin + 0.5*DY) + YINDEX(NUMROWS)));
[X,Y] = meshgrid(xloc,yloc);
[xX,xY] = meshgrid(x_xloc,yloc);
[yX,yY] = meshgrid(xloc,y_yloc);

[u1,v2,u3,v4] = aLLOCfACE(u,v);
UAv = 0.5.*(u1 + u3);
VAv = 0.5.*(v2 + v4);
NVel = sqrt(UAv.^2 + VAv.^2);
for i=1:NUMROWS
   for j=1:NUMCOLS
      CalVel(i,j) = NVel(((i-1)*NUMCOLS)+j);
      CalDep(i,j) = H(((i-1)*NUMCOLS)+j);
      PVv(i,j) = v(((i-1)*NUMCOLS)+j);
      fprintf(f1,'%10.3f\t %11.3f\t %8.4f\n',xloc(j),yloc(i),CalVel(i,j));
      fprintf(f2,'%10.3f\t %11.3f\t %8.4f\n',xloc(j),yloc(i),CalDep(i,j));
   end
end
PVv(NUMROWS+1,:) = v(ROWEND(NUMROWS)+1:1:NUMINCY)';
for i=1:NUMROWS
   for j=1:NUMCOLS+1
      PUv(i,j) = u(((i-1)*XINC)+j);
   end
end

fclose(f1);
fclose(f2);

d5 = 0.0:2.0:20.0;
v5 = [0.1 0.3 0.5 0.7 0.8 0.9 1.0 1.1 1.2 1.3];
v8 = [-2.0 -1.5 -1.0 -0.5 -0.3 -0.1 0.1 0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7 1.9 2.1];

h2 = figure(2);
% top left.
subplot(2,2,1);
hold on;
[D1,k1] = contour(X,Y,CalVel,v5,'-k');
e1 = clabel(D1);
title('(a)','FontSize',12);
xlabel('Easting in meters','FontSize',12);
ylabel('Northing in meters','FontSize',12);
axis([XMin XMax YMin YMax]);
pbaspect([2 1 1]);
set(gca,'XTick',PlotX);
set(gca,'YTick',PlotY);
set(gca,'XTickLabel',int2str(PlotX'),'FontSize',8);
set(gca,'YTickLabel',int2str(PlotY'),'FontSize',8);
set(gca,'XGrid','on','GridLineStyle',':');
set(gca,'YGrid','on','GridLineStyle',':');
hold off;
% top right.
subplot(2,2,2);
hold on;
[D2,k2] = contour(X,Y,CalDep,d5,'-k');
e2 = clabel(D2);
title('(b)','FontSize',12);
xlabel('Easting in meters','FontSize',12);
ylabel('Northing in meters','FontSize',12);
axis([XMin XMax YMin YMax]);
pbaspect([2 1 1]);
set(gca,'XTick',PlotX);
set(gca,'YTick',PlotY);
set(gca,'XTickLabel',int2str(PlotX'),'FontSize',8);
set(gca,'YTickLabel',int2str(PlotY'),'FontSize',8);
set(gca,'XGrid','on','GridLineStyle',':');
set(gca,'YGrid','on','GridLineStyle',':');
hold off;
% bottom left.
subplot(2,2,3);
hold on;
[D3,k3] = contour(xX,xY,PUv,v8,'-k');
e3 = clabel(D3);
title('(c)','FontSize',12);
xlabel('Easting in meters','FontSize',12);
ylabel('Northing in meters','FontSize',12);
axis([XMin XMax YMin YMax]);
pbaspect([2 1 1]);
set(gca,'XTick',PlotX);
set(gca,'YTick',PlotY);
set(gca,'XTickLabel',int2str(PlotX'),'FontSize',8);
set(gca,'YTickLabel',int2str(PlotY'),'FontSize',8);
set(gca,'XGrid','on','GridLineStyle',':');
set(gca,'YGrid','on','GridLineStyle',':');
hold off;
% bottom right.
subplot(2,2,4);
hold on;
[D4,k4] = contour(yX,yY,PVv,v8,'-k');
e4 = clabel(D4);
title('(d)','FontSize',12);
xlabel('Easting in meters','FontSize',12);
ylabel('Northing in meters','FontSize',12);
axis([XMin XMax YMin YMax]);
pbaspect([2 1 1]);
set(gca,'XTick',PlotX);
set(gca,'YTick',PlotY);
set(gca,'XTickLabel',int2str(PlotX'),'FontSize',8);
set(gca,'YTickLabel',int2str(PlotY'),'FontSize',8);
set(gca,'XGrid','on','GridLineStyle',':');
set(gca,'YGrid','on','GridLineStyle',':');
hold off;
pause(2);
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.50 1.5 10.0 6.5]);
saveas(gcf,fpeps,'ps2');
saveas(gcf,fpmat,'fig');

clear D CalDepth CalVel e NVel UAv v8 v5 VAv X xloc Y yloc PUvj;
clear PVv x_xloc y_yloc;
clear xX xY yX yY d5 fnv fnd f1 f2 PlotX PlotY fpeps fpmat;
clear XMin XMax YMin YMax;
clear D1 D2 D3 D4 k1 k2 k3 k4 e1 e2 e3 e4;
clear u1 u3 v2 v4;
close(h2);
return;
%EOF
