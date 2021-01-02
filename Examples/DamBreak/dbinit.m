function dbinit
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% dbinit sets things up for tracking of water depth at prescribed
% locations for simulations of the dambreak style flume experiment of 
% Bellos et al. (1992).

global DepthLoc1 DepthLoc2 DepthLoc3 DepthLoc4 DepthLoc5 DepthLoc6 DepthLoc7
global DepthLoc8 ENDTIME FLUID_DT H STARTTIME TimeCounter
global L1Size L2Size L3Size L4Size L5Size L6Size L7Size L8Size L1Vol L2Vol
global L3Vol L4Vol L5Vol L6Vol L7Vol L8Vol 

%
TotTime = double(0.0);
NumInc = 0;
%
TotTime = (ENDTIME - STARTTIME)*(24.0*60.0*60.0);
NumInc = floor(TotTime/FLUID_DT);

% initialize storage vectors.
DepthLoc1 = zeros(NumInc,1);
DepthLoc2 = zeros(NumInc,1);
DepthLoc3 = zeros(NumInc,1);
DepthLoc4 = zeros(NumInc,1);
DepthLoc5 = zeros(NumInc,1);
DepthLoc6 = zeros(NumInc,1);
DepthLoc7 = zeros(NumInc,1);
DepthLoc8 = zeros(NumInc,1);
TimeCounter = zeros(NumInc,1);

L1Size = 2;
L1Vol = [ 45 46 ]; 
L2Size = 4;
L2Vol = [ 1095 1096 1125 1126 ];
L3Size = 2;
L3Vol = [ 2047 2048 ];
L4Size = 2;
L4Vol = [ 2077 2078 ];
L5Size = 2;
L5Vol = [ 2649 2679 ];
L6Size = 2;
L6Vol = [ 3253 3283 ];
L7Size = 4;
L7Vol = [ 3855 3856 3885 3886 ];
L8Size = 4;
L8Vol = [ 4455 4456 4485 4486 ];

% Calcs.
if (L1Size > 0)
   Dep = double(0.0);
   for i=1:L1Size
       Dep = Dep + H(L1Vol(i));
   end
   Dep = Dep*(1/L1Size);
   DepthLoc1(1,1) = Dep;
end
if (L2Size > 0)
   Dep = double(0.0);
   for i=1:L2Size
       Dep = Dep + H(L2Vol(i));
   end
   Dep = Dep*(1/L2Size);
   DepthLoc2(1,1) = Dep;
end
if (L3Size > 0)
   Dep = double(0.0);
   for i=1:L3Size
       Dep = Dep + H(L3Vol(i));
   end
   Dep = Dep*(1/L3Size);
   DepthLoc3(1,1) = Dep;
end
if (L4Size > 0)
   Dep = double(0.0);
   for i=1:L4Size
       Dep = Dep + H(L4Vol(i));
   end
   Dep = Dep*(1/L4Size);
   DepthLoc4(1,1) = Dep;
end
if (L5Size > 0)
   Dep = double(0.0);
   for i=1:L5Size
       Dep = Dep + H(L5Vol(i));
   end
   Dep = Dep*(1/L5Size);
   DepthLoc5(1,1) = Dep;
end
if (L6Size > 0)
   Dep = double(0.0);
   for i=1:L6Size
       Dep = Dep + H(L6Vol(i));
   end
   Dep = Dep*(1/L6Size);
   DepthLoc6(1,1) = Dep;
end
if (L7Size > 0)
   Dep = double(0.0);
   for i=1:L7Size
       Dep = Dep + H(L7Vol(i));
   end
   Dep = Dep*(1/L7Size);
   DepthLoc7(1,1) = Dep;
end
if (L8Size > 0)
   Dep = double(0.0);
   for i=1:L8Size
       Dep = Dep + H(L8Vol(i));
   end
   Dep = Dep*(1/L8Size);
   DepthLoc8(1,1) = Dep;
end

clear Dep i NumInc TotTime;
return;
%EOF
