function [A,B,C,D] = ss(Hd)
%SS  Discrete-time filter to state-space conversion.
%   [A,B,C,D] = SS(Hd) converts discrete-time filter Hd to state-space
%   representation given by 
%     x(k+1) = A*x(k) + B*u(k)
%     y(k)   = C*x(k) + D*u(k)
%   where x is the state vector, u is the input vector, and y is the output
%   vector. 
%
%   See also DFILT.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.7 $  $Date: 2002/05/04 01:48:56 $

k = Hd.Lattice;

% Find the state-space of the min-phase, we use a and b only
[A,B] = super_ss(Hd);

% Make sure k is a row
k = k(:).';

% Force uniformity with empty state matrix.
if isempty(A)
  A = [];
  B = zeros(0,1);
  C = zeros(1,0);
else
  C = [k(end)'.*k(1:end-1) 1];
end

D = k(end)';

