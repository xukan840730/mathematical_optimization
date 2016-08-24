function [SD, dirType] = compdir(Z,H,gf,nvars,f)
% COMPDIR Computes a search direction in a subspace defined by Z. 
%    Helper function for NLCONST.
%    Returns Newton direction if possible.
%    Returns random direction if gradient is small.
%    Otherwise, returns steepest descent direction.  
%    If the steepest descent direction is small it computes a negative
%    curvature direction based on the most negative eigenvalue.
%    For singular matrices, returns steepest descent even if small.

%   Copyright 1990-2002 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2002/03/29 18:32:00 $
%   Mary Ann Branch 10-20-96.

% Define constant strings
Newton = 'Newton';
Random = 'random';
SteepDescent = 'steepest descent';
Eigenvector = 'eigenvector';

%  SD=-Z*((Z'*H*Z)\(Z'*gf));
dirType = [];
% Compute the projected Newton direction if possible
 projH = Z'*H*Z;
 [R, p] = chol(projH);
 if ~p  % positive definite: use Newton direction
   SD = - Z*(R \ ( R'\(Z'*gf)));
   dirType = Newton;
 else % not positive definite
   % If the gradient is small, try a random direction:
      % Sometimes the search direction goes to zero in negative
      % definite problems when the current point rests on
      % the top of the quadratic function. In this case we can move in
      % any direction to get an improvement in the function so 
      % foil search direction by giving a random gradient.
   if norm(gf) < sqrt(eps)
      SD = -Z*Z'*(rand(nvars,1) - 0.5);
      dirType = Random;
   else
     % steepest descent
     stpDesc = - Z*(Z'*gf);
     % check if ||SD|| is close to zero 
     if norm(stpDesc) > sqrt(eps)       
        SD = stpDesc;
        dirType = SteepDescent;
     else
        % Look for a negative curvature direction
        %  Some attempt at efficiency: usually it's
        %  faster to use EIG unless many variables.
        if nvars < 400  
           [VV,DD] = eig(projH);
           [smallRealEig, eigind] = min(diag(DD));
           ev = VV(:,eigind(1));
        else
           options.disp = 0;
           [ev, smallRealEig, flag] = eigs(projH,1,'sr',options);
           if flag  % Call to eigs failed
              [VV,DD] = eig(projH);
              [smallRealEig, eigind] = min(diag(DD));
              ev = VV(:,eigind(1));
           end
        end
        
        if smallRealEig < 0
          % check the sign of SD and the magnitude.
          SDtol = 100*eps*norm(gf); % Note: we know norm(gf) > sqrt(eps)
          Zev = Z*ev;
          if Zev'*gf > SDtol
            SD = -Zev;
            dirType = Eigenvector;
          elseif Zev'*gf < SDtol
            SD = Zev;
            dirType = Eigenvector;
          else % 
            SD = stpDesc;
            dirType = SteepDescent; 
          end
        else % The projected Hessian is singular,i.e., zero direction is ok
          %  -- will propagate thru the algorithm.
          SD = stpDesc;
          dirType = SteepDescent; 
        end % smallRealEig < 0
      end % randSD'*(gf) < -SDtol
   end %  norm(stpDesc) > sqrt(eps)  
 end % ~p

   
