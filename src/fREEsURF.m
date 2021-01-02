function FSurfNew = fREEsURF(FSurfNew,NowTime,ocntr,cETime)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % fREEsURF calculates the new free surface elevations across the entire domain.
    % This function will set-up both the right-hand side (RHS) and the left-hand
    % side (LHS) of the system of equations for free surface elevations.  Free
    % surface elevations are solved in a semi-implicit manner using pre-conditioned
    % conjugate gradients.
    %
    % Received:
    %
    % FSurfNew = EtaNew [NUMNODES,1] new free surface elevations.
    % NowTime [1] = current time in simulation. (fluid_t).
    % ocntr [1] = output counter value.
    % cETime (float) = current elapsed time for time series lookup
    % 
    % Returned:
    % FSurfNew = EtaNew [NUMNODES,1] new free surface elevations.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Copyright and License
    %
    % Copyright 2021 Nick Martin
    %
    % This file is part of MOD_FreeSurf2D.
    %
    % MOD_FreeSurf2d is free software: you can redistribute it and/or modify
    % it under the terms of the GNU Affero General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    %
    % MOD_FreeSurf2D is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU Affero General Public License for more details.
    %
    % You should have received a copy of the GNU Affero General Public License
    % along with MOD_FreeSurf2D.  If not, see <https://www.gnu.org/licenses/>.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global EPSILON MAXITER NUMNODES OUTINT OUTP PRECOND qRHS sLHS

    % local variables.
    Delta = double(zeros(NUMNODES,1));      % Delta vector for RHS calculation.
    EToEta = double(zeros(NUMNODES,1));     % Normalizing tranformation factor.
    SFlag = 0;                              % Flag for pcg -- integer.
    SRelRes = double(0.0);                  % pcg residual.
    SIter = 0;                              % Iterations at convergence.
    x0 = double(ones(NUMNODES,1));          % potential initial guess for pcg.

    %Now calculate the values for the RHS of the free surface solution.
    Delta = dELTAcALC(Delta);
    % Get the RHS.
    qRHS = rhscALC(qRHS,Delta);
    % Now setup the LHS.
    % return EToEta and qRHS but also return sLHS as a global for memory
    %  conservation.
    [EToEta,qRHS] = lhscALC(EToEta,qRHS,cETime);
    % Conjugate gradient solver.
    if (PRECOND == 0)
       [FSurfNew,SFlag,SRelRes,SIter] = pcg(sLHS,qRHS,EPSILON,MAXITER,[],[],x0);
    elseif (PRECOND == 1)
       [FSurfNew,SFlag,SRelRes,SIter] = pcg(sLHS,qRHS,EPSILON,MAXITER,...
          spdiags(diag(sLHS,0),0,NUMNODES,NUMNODES),[],x0);   
    elseif (PRECOND == 2)
       R = ichol(sLHS, struct('type','ict','droptol',1e-3,'shape','upper'));
       [FSurfNew,SFlag,SRelRes,SIter] = pcg(sLHS,qRHS,EPSILON,MAXITER,R,R',x0);
       clear R;
    else
       fprintf(OUTP,'Invalid value for PRECOND in mainfluid.m \n');
       disp('Invalid value for PRECOND in mainfluid.m');
    end
    % Adjust values from normalized.
    FSurfNew = FSurfNew.*EToEta;
    % Output to file
    if (ocntr == OUTINT)
       fprintf('Current time = %20.8f [days]\n',cETime);
       fprintf('Flag = %1i\t Iter. = %4i\t Resid.= %12.3E\n',...
               SFlag,SIter,SRelRes);
    end
    % not needed
    %clear Delta EToEta SFlag SRelRes SIter x0;
    %clear NowTime ocntr;
end

%EOF