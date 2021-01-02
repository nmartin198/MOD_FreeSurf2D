function dbset(Counter,PTime)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% dbset saves water depths and time steps for the dam break test.  The depth
% averaged to the center of each volume (H) is employed for the depth 
% calculations.  Those locations at a volume corner are obtained by 
% bilinear interpolation. This script is for simulation of the dambreak style
% flume experiment of Bellos et al. (1992).

global DepthLoc1 DepthLoc2 DepthLoc3 DepthLoc4 DepthLoc5 DepthLoc6 DepthLoc7
global DepthLoc8 H TimeCounter
global L1Size L2Size L3Size L4Size L5Size L6Size L7Size L8Size L1Vol L2Vol
global L3Vol L4Vol L5Vol L6Vol L7Vol L8Vol 

% Initialize.
NI = size(TimeCounter,1);
% Calculations.
if (Counter <= NI)
   TimeCounter(Counter,1) = PTime;
   if (L1Size > 0)
      Dep = double(0.0);
      for i=1:L1Size
         Dep = Dep + H(L1Vol(i));
      end
      Dep = Dep*(1/L1Size);
      DepthLoc1(Counter,1) = Dep;
   end
   if (L2Size > 0)
      Dep = double(0.0);
      for i=1:L2Size
         Dep = Dep + H(L2Vol(i));
      end
      Dep = Dep*(1/L2Size);
      DepthLoc2(Counter,1) = Dep;
   end
   if (L3Size > 0)
      Dep = double(0.0);
      for i=1:L3Size
         Dep = Dep + H(L3Vol(i));
      end
      Dep = Dep*(1/L3Size);
      DepthLoc3(Counter,1) = Dep;
   end
   if (L4Size > 0)
      Dep = double(0.0);
      for i=1:L4Size
         Dep = Dep + H(L4Vol(i));
      end
      Dep = Dep*(1/L4Size);
      DepthLoc4(Counter,1) = Dep;
   end
   if (L5Size > 0)
      Dep = double(0.0);
      for i=1:L5Size
         Dep = Dep + H(L5Vol(i));
      end
      Dep = Dep*(1/L5Size);
      DepthLoc5(Counter,1) = Dep;
   end
   if (L6Size > 0)
      Dep = double(0.0);
      for i=1:L6Size
         Dep = Dep + H(L6Vol(i));
      end
      Dep = Dep*(1/L6Size);
      DepthLoc6(Counter,1) = Dep;
   end
   if (L7Size > 0)
      Dep = double(0.0);
      for i=1:L7Size
         Dep = Dep + H(L7Vol(i));
      end
      Dep = Dep*(1/L7Size);
      DepthLoc7(Counter,1) = Dep;
   end
   if (L8Size > 0)
      Dep = double(0.0);
      for i=1:L8Size
         Dep = Dep + H(L8Vol(i));
      end
      Dep = Dep*(1/L8Size);
      DepthLoc8(Counter,1) = Dep;
   end
end

clear Dep i Counter PTime NI;
return;
%EOF
