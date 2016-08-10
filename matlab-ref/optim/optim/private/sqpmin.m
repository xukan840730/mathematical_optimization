function [x,val,output,exitflag,lambda] = sqpmin(c,H,mtxmpy,xstart,A,b,lb,ub,verb,...
    options,defaultopt,computeLambda,computeConstrViolation,varargin)
%SQPMIN	Solve quadratic problems with box constraints or linear equalities
%
% Locate local solution to
%
%        min { q(x) = .5x'Hx + c'x : lb <= x <= ub}. 
%
%                            or
%
%       min { q(x) = .5x'Hx + c'x : Ax = b},
%
%
% where H is a sparse symmetric matrix. 

%   Copyright 1990-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2009/05/07 18:25:25 $

if nargin < 2
  error('optim:sqpmin:NotEnoughInputs','sqpmin requires at least 2 arguments.')
end
if nargin <=2, xstart = []; end
n = length(c); 

if isempty(lb), lb = -inf*ones(n,1); end
if isempty(ub),ub = inf*ones(n,1); end
arg = (ub >= 1e10); arg2 = (lb <= -1e10);
ub(arg) = inf;
lb(arg2) = -inf;
if any(ub == lb) 
   error('optim:sqpmin:EqualLowerUpperBnd', ...
         ['Equal upper and lower bounds not permitted in this large-scale algorithm.\n' ...
          'Use equality constraints and the medium-scale algorithm instead.'])
elseif min(ub-lb) <= 0
   error('optim:sqpmin:InconsistentBnds','Inconsistent bounds.')
end
if isempty(xstart) 
    xstart = startx(ub,lb); 
end
% If components of initial x not within bounds, set those components
% of initial point to a "box-centered" point
xinitOutOfBounds_idx = xstart < lb | xstart > ub;
if any(xinitOutOfBounds_idx)
    xstart = startx(ub,lb,xstart,xinitOutOfBounds_idx);
end

numberOfVariables = n;
if n == 0
  error('optim:sqpmin:InvalidN','n must be positive.')
end

%   INITIALIZATIONS
lambda.lower = [];
lambda.upper = [];
lambda.eqlin = [];  
lambda.ineqlin = [];  % This won't change because no inequalities.
val = []; gopt=[];
output = [];
it = 1; 
if nargin < 5, A = []; end
if isempty(A)
   % Box-constrained problem
   % Set default values for TolFun and MaxIter
   defaultopt.TolFun = 100*eps;
   defaultopt.MaxIter = 200;   

   [x,val,gopt,it,npcg,exitflag,lambda,msg] = sqpbox(c,H,mtxmpy,lb,ub,xstart, ...
       options,defaultopt,numberOfVariables,verb,computeLambda,varargin{:});
  lambda.ineqlin = []; lambda.eqlin = [];
  output.iterations = it;   
  if computeConstrViolation
      output.constrviolation = max([0; (lb-x); (x-ub) ]);
  else
      output.constrviolation = [];
  end
  output.algorithm = 'large-scale: reflective trust-region';
  output.firstorderopt = gopt;
  output.cgiterations = npcg;
  output.message = msg;
else
   if (max(lb) > -inf) || (min(ub) < inf)
      error('optim:sqpmin:InvalidConstraints', ...
            'sqpmin doesn''t handle both box constraints and Ax = b.');
   else
      % Equality constrained problem
      [mA,nA] = size(A);
      if nargin < 6, b = zeros(mA,1); end
      if isempty(b), b = zeros(mA,1); end
      
      % Since we only call PPCGR once, use the original problem tolerances,
      % TolFun and MaxIter, instead of the inner iteration tolerances,
      % TolPCG and MaxPCGIter.
            
      % Set default values for TolFun and MaxIter
      defaultopt.TolFun = 1e-6;
      defaultopt.MaxIter = 2*max((numberOfVariables - mA),1);
      tolPCG     = optimget(options,'TolPCG',defaultopt,'fast');
      maxPCGiter = optimget(options,'MaxPCGIter',defaultopt,'fast');
      tolFun     = optimget(options,'TolFun',defaultopt,'fast');
      maxIter    = optimget(options,'MaxIter',defaultopt,'fast');
      
      if tolPCG ~= defaultopt.TolPCG
         warning('optim:sqpmin:TolPCGignored', ...
             ['TolPCG is ignored for large-scale, equality constrained ', ...
             'problems. Use TolFun instead.']);
      end
      if maxPCGiter ~= defaultopt.MaxPCGIter
         warning('optim:sqpmin:MaxPCGIterignored', ...
             ['MaxPCGIter is ignored for large-scale, equality constrained ', ...
             'problems. Use MaxIter instead.']);
      end

      % Here we must be careful not to create an options structure if one does not
      % already exist. This will lead to an error in ppcgr.m when trying to read
      % options.Preconditioner, an internal option. In this case, just update defaultopt.
      if isempty(options)
         defaultopt.TolPCG = tolFun;
         defaultopt.MaxPCGIter = maxIter;
      else
         options.TolPCG = tolFun;
         options.MaxPCGIter = maxIter;
      end
      
      % Note we pass options in so some values are different for PPCGR than for SQPBOX
      [x,po,npcg,pgnrm,exitflag,lambda,msg]=ppcgr(c,H,mtxmpy,A,b,options,defaultopt,verb,computeLambda,...,
          [],[],[],[],[],[],varargin{:});
      if exitflag == -10
         % ppcgr aborted
         return
      end
      
      it = 1; % number of iterations reported in output.cgiterations
      % Calculate final objective value 
      w = feval(mtxmpy,H,x,varargin{:}); 
      gopt = pgnrm;
      val = x'*(c + .5*w);
      
      output.iterations = it; 
      if computeConstrViolation
          output.constrviolation = norm(A*x-b, inf);
      else
          output.constrviolation = [];
      end
      output.algorithm = 'large-scale: projective preconditioned conjugate gradients';
      output.firstorderopt = gopt;
      output.cgiterations = npcg;
      output.message = msg;
   end
end





