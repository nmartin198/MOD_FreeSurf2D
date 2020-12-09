function rEADiNPUT
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% rEADiNPUT reads the parameters from the input file input.txt.

global CENTRALLATITUDE COROMEGA DATUM FLUID_DT DX DY ENDTIME
global EPSILON EVIS G GAMMATX GAMMATY HCUTOFF KAPPA MAXCR MAXITER
global MAXSTEPS MINSTEPS MNWAL NUK NUMCOLS NUMROWS OUTINT 
global PATHTRAC PREC PRECH PRECOND QINBC QINXFLUX QINXVOL QINYFLUX
global QINYVOL RADFLUXBC RADORLFSBC RADVELBC RFLUXXFLUX RFLUXYFLUX
global RORLFSXVOL RORLFSYVOL RVELXVOL RVELYVOL RHOW SIMNAME
global STARTTIME TDEPDIRCBC TDEPDXDEP TDEPDXVOL TDEPDYDEP TDEPDYVOL
global THETA UA VA VELDIRCBC VELDXVOL VELDYVOL VELDXVEL
global VELDYVEL

% Initialize all of the non-string globals and all of the non-source
%  values.
PREC = double(0.0);
PRECH = double(0.0);
STARTTIME = double(0.0);
ENDTIME = double(0.0);
DATUM = double(0.0);
DX = double(0.0);
DY = double(0.0);
NUMROWS = 0;
NUMCOLS = 0;
OUTINT = 0;
FLUID_DT = double(0.0);
THETA = double(0.0);
HCUTOFF = double(0.0);
EPSILON = double(0.0);
MAXITER = 0;
PRECOND = 0;
MAXCR = 0;
PATHTRAC = 0;
MAXSTEPS = 0;
MINSTEPS = 0;
G = double(0.0);
KAPPA = double(0.0);
MNWAL = double(0.0);
NUK = double(0.0);
RHOW = double(0.0);
EVIS = double(0.0);
CENTRALLATITUDE = double(0.0);
COROMEGA = double(0.0);
GAMMATX = double(0.0);
UA = double(0.0);
GAMMATY = double(0.0);
VA = double(0.0);
TDEPDIRCBC = 0;
VELDIRCBC = 0;
QINBC = 0;
RADVELBC = 0;
RADORLFSBC = 0;
RADFLUXBC = 0;
% local variable.
INP = 0;

% open input file.
INP = fopen('input.txt');
% read through file looking for parameters.
Liner = fgetl(INP);
while (Liner(1) ~= -1)
   if (Liner(1) ~= '#') 
      if (~(isempty(findstr(Liner,'STARTTIME ='))))
         STARTTIME = sscanf(Liner,'STARTTIME =%f',1);
      end
      if (~(isempty(findstr(Liner,'ENDTIME ='))))
         ENDTIME = sscanf(Liner,'ENDTIME =%f',1);
      end
      if (~(isempty(findstr(Liner,'SIMNAME ='))))
         SIMNAME = sscanf(Liner,'SIMNAME =%s',1);
      end
      if (~(isempty(findstr(Liner,'DATUM ='))))
         DATUM = sscanf(Liner,'DATUM =%f',1);
      end
      if (~(isempty(findstr(Liner,'DX ='))))
         DX = sscanf(Liner,'DX =%f',1);
      end
      if (~(isempty(findstr(Liner,'DY ='))))
         DY = sscanf(Liner,'DY =%f',1);
      end
      if (~(isempty(findstr(Liner,'NUMROWS ='))))
         NUMROWS = sscanf(Liner,'NUMROWS =%i',1);
      end
      if (~(isempty(findstr(Liner,'NUMCOLS ='))))
         NUMCOLS = sscanf(Liner,'NUMCOLS =%i',1);
      end
      if (~(isempty(findstr(Liner,'OUTINT ='))))
         OUTINT = sscanf(Liner,'OUTINT =%i',1);
      end
      if (~(isempty(findstr(Liner,'FLUID_DT ='))))
         FLUID_DT = sscanf(Liner,'FLUID_DT =%f',1);
      end
      if (~(isempty(findstr(Liner,'THETA ='))))
         THETA = sscanf(Liner,'THETA =%f',1);
      end
      if (~(isempty(findstr(Liner,'HCUTOFF ='))))
         HCUTOFF = sscanf(Liner,'HCUTOFF =%f',1);
      end
      if (~(isempty(findstr(Liner,'EPSILON ='))))
         EPSILON = sscanf(Liner,'EPSILON =%E',1);
      end
      if (~(isempty(findstr(Liner,'MAXITER ='))))
         MAXITER = sscanf(Liner,'MAXITER =%i',1);
      end
      if (~(isempty(findstr(Liner,'PRECOND ='))))
         PRECOND = sscanf(Liner,'PRECOND =%i',1);
      end
      if (~(isempty(findstr(Liner,'MAXCR ='))))
         MAXCR = sscanf(Liner,'MAXCR =%i',1);
      end
      if (~(isempty(findstr(Liner,'PATHTRAC ='))))
         PATHTRAC = sscanf(Liner,'PATHTRAC =%i',1);
      end
      if (~(isempty(findstr(Liner,'MAXSTEPS ='))))
         MAXSTEPS = sscanf(Liner,'MAXSTEPS =%i',1);
      end
      if (~(isempty(findstr(Liner,'MINSTEPS ='))))
         MINSTEPS = sscanf(Liner,'MINSTEPS =%i',1);
      end
      if (~(isempty(findstr(Liner,'G ='))))
         G = sscanf(Liner,'G =%f',1);
      end
      if (~(isempty(findstr(Liner,'KAPPA ='))))
         KAPPA = sscanf(Liner,'KAPPA =%f',1);
      end
      if (~(isempty(findstr(Liner,'MNWAL ='))))
         MNWAL = sscanf(Liner,'MNWAL =%f',1);
      end
      if (~(isempty(findstr(Liner,'NUK ='))))
         NUK = sscanf(Liner,'NUK =%E',1);
      end
      if (~(isempty(findstr(Liner,'RHOW ='))))
         RHOW = sscanf(Liner,'RHOW =%f',1);
      end
      if (~(isempty(findstr(Liner,'EVIS ='))))
         EVIS = sscanf(Liner,'EVIS =%E',1);
      end
      if (~(isempty(findstr(Liner,'CENTRALLATITUDE ='))))
         CENTRALLATITUDE = sscanf(Liner,'CENTRALLATITUDE =%f',1);
      end
      if (~(isempty(findstr(Liner,'COROMEGA ='))))
         COROMEGA = sscanf(Liner,'COROMEGA =%f',1);
      end
      if (~(isempty(findstr(Liner,'GAMMATX ='))))
         GAMMATX = sscanf(Liner,'GAMMATX =%f',1);
      end
      if (~(isempty(findstr(Liner,'UA ='))))
         UA = sscanf(Liner,'UA =%f',1);
      end
      if (~(isempty(findstr(Liner,'GAMMATY ='))))
         GAMMATY = sscanf(Liner,'GAMMATY =%f',1);
      end
      if (~(isempty(findstr(Liner,'VA ='))))
         VA = sscanf(Liner,'VA =%f',1);
      end
      if (~(isempty(findstr(Liner,'TDEPDIRCBC ='))))
         TDEPDIRCBC = sscanf(Liner,'TDEPDIRCBC =%i',1);
      end
      if (~(isempty(findstr(Liner,'TDEPDXVOL ='))))
         if (TDEPDIRCBC ~= 0)
            be = findstr(Liner,'[');
            en = findstr(Liner,']');
            TDEPDXVOL = (sscanf(Liner(be+1:en-1),'%i',inf))';
            clear be en;
         else
            TDEPDXVOL = 0;
         end
      end
      if (~(isempty(findstr(Liner,'TDEPDXDEP ='))))
         if (TDEPDIRCBC ~= 0)
            sizer = size(TDEPDXVOL,2);
            be = findstr(Liner,'=');
            en = size(Liner,2);
            TDEPDXDEP = (sscanf(Liner(be+1:en),'%f',sizer))';
            clear be en sizer;
         else
            TDEPDXDEP = 0.0;
         end
      end
      if (~(isempty(findstr(Liner,'TDEPDYVOL ='))))
         if (TDEPDIRCBC ~= 0)
            be = findstr(Liner,'[');
            en = findstr(Liner,']');
            TDEPDYVOL = (sscanf(Liner(be+1:en-1),'%i',inf))';
            clear be en;
         else
            TDEPDYVOL = 0;
         end
      end
      if (~(isempty(findstr(Liner,'TDEPDYDEP ='))))
         if (TDEPDIRCBC ~= 0)
            sizer = size(TDEPDYVOL,2);
            be = findstr(Liner,'=');
            en = size(Liner,2);
            TDEPDYDEP = (sscanf(Liner(be+1:en),'%f',sizer))';
            clear be en sizer;
         else
            TDEPDYDEP = 0.0;
         end
      end
      if (~(isempty(findstr(Liner,'VELDIRCBC ='))))
         VELDIRCBC = sscanf(Liner,'VELDIRCBC =%i',1);
      end
      if (~(isempty(findstr(Liner,'VELDXVOL ='))))
         if (VELDIRCBC ~= 0)
            be = findstr(Liner,'[');
            en = findstr(Liner,']');
            VELDXVOL = (sscanf(Liner(be+1:en-1),'%i',inf))';
            clear be en;
         else
            VELDXVOL = 0;
         end
      end
      if (~(isempty(findstr(Liner,'VELDXVEL ='))))
         if (VELDIRCBC ~= 0)
            sizer = size(VELDXVOL,2);
            be = findstr(Liner,'=');
            en = size(Liner,2);
            VELDXVEL = (sscanf(Liner(be+1:en),'%f',sizer))';
            clear be en sizer;
         else
            VELDXVEL = 0.0;
         end
      end
      if (~(isempty(findstr(Liner,'VELDYVOL ='))))
         if (VELDIRCBC ~= 0)
            be = findstr(Liner,'[');
            en = findstr(Liner,']');
            VELDYVOL = (sscanf(Liner(be+1:en-1),'%i',inf))';
            clear be en;
         else
            VELDYVOL = 0;
         end
      end
      if (~(isempty(findstr(Liner,'VELDYVEL ='))))
         if (VELDIRCBC ~= 0)
            sizer = size(VELDYVOL,2);
            be = findstr(Liner,'=');
            en = size(Liner,2);
            VELDYVEL = (sscanf(Liner(be+1:en),'%f',sizer))';
            clear be en sizer;
         else
            VELDYVEL = 0.0;
         end
      end
      if (~(isempty(findstr(Liner,'QINBC ='))))
         QINBC = sscanf(Liner,'QINBC =%i',1);
      end
      if (~(isempty(findstr(Liner,'QINXVOL ='))))
         if (QINBC ~= 0)
            be = findstr(Liner,'[');
            en = findstr(Liner,']');
            QINXVOL = (sscanf(Liner(be+1:en-1),'%i',inf))';
            clear be en;
         else
            QINXVOL = 0;
         end
      end
      if (~(isempty(findstr(Liner,'QINXFLUX ='))))
         if (QINBC ~= 0)
            sizer = size(QINXVOL,2);
            be = findstr(Liner,'=');
            en = size(Liner,2);
            QINXFLUX = (sscanf(Liner(be+1:en),'%f',sizer))';
            clear be en sizer;
         else
            QINXFLUX = 0.0;
         end
      end
      if (~(isempty(findstr(Liner,'QINYVOL ='))))
         if (QINBC ~= 0)
            be = findstr(Liner,'[');
            en = findstr(Liner,']');
            QINYVOL = (sscanf(Liner(be+1:en-1),'%i',inf))';
            clear be en;
         else
            QINYVOL = 0;
         end
      end
      if (~(isempty(findstr(Liner,'QINYFLUX ='))))
         if (QINBC ~= 0)
            sizer = size(QINYVOL,2);
            be = findstr(Liner,'=');
            en = size(Liner,2);
            QINYFLUX = (sscanf(Liner(be+1:en),'%f',sizer))';
            clear be en sizer;
         else
            QINYFLUX = 0.0;
         end
      end
      if (~(isempty(findstr(Liner,'RADVELBC ='))))
         RADVELBC = sscanf(Liner,'RADVELBC =%i',1);
      end
      if (~(isempty(findstr(Liner,'RVELXVOL ='))))
         if (RADVELBC ~= 0)
            be = findstr(Liner,'[');
            en = findstr(Liner,']');
            RVELXVOL = (sscanf(Liner(be+1:en-1),'%i',inf))';
            clear be en;
         else
            RVELXVOL = 0;
         end
      end
      if (~(isempty(findstr(Liner,'RVELYVOL ='))))
         if (RADVELBC ~= 0)
            be = findstr(Liner,'[');
            en = findstr(Liner,']');
            RVELYVOL = (sscanf(Liner(be+1:en-1),'%i',inf))';
            clear be en;
         else
            RVELYVOL = 0;
         end
      end
      if (~(isempty(findstr(Liner,'RADORLFSBC ='))))
         RADORLFSBC = sscanf(Liner,'RADORLFSBC =%i',1);
      end
      if (~(isempty(findstr(Liner,'RORLFSXVOL ='))))
         if (RADORLFSBC ~= 0)
            be = findstr(Liner,'[');
            en = findstr(Liner,']');
            RORLFSXVOL = (sscanf(Liner(be+1:en-1),'%i',inf))';
            clear be en;
         else
            RORLFSXVOL = 0;
         end
      end
      if (~(isempty(findstr(Liner,'RORLFSYVOL ='))))
         if (RADORLFSBC ~= 0)
            be = findstr(Liner,'[');
            en = findstr(Liner,']');
            RORLFSYVOL = (sscanf(Liner(be+1:en-1),'%i',inf))';
            clear be en;
         else
            RORLFSYVOL = 0;
         end
      end
      if (~(isempty(findstr(Liner,'RADFLUXBC ='))))
         RADFLUXBC = sscanf(Liner,'RADFLUXBC =%i',1);
      end
      if (~(isempty(findstr(Liner,'RFLUXXFLUX ='))))
         if (RADFLUXBC ~= 0)
            RFLUXXFLUX = (sscanf(Liner,'RFLUXXFLUX =%f',1))';
            clear be en;
         else
            RFLUXXFLUX = 0.0;
         end
      end
      if (~(isempty(findstr(Liner,'RFLUXYFLUX ='))))
         if (RADFLUXBC ~= 0)
            RFLUXYFLUX = (sscanf(Liner,'RFLUXYFLUX =%f',1))';
            clear be en;
         else
            RFLUXYFLUX = 0.0;
         end
      end
   end
   Liner = fgetl(INP);
end

fclose(INP);
clear Liner INP;
return;
%EOF