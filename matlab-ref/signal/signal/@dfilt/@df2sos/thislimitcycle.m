function [y,zi,overflows] = thislimitcycle(Hd,x)
%THISLIMITCYCLE   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2005/10/14 16:24:58 $

q = Hd.filterquantizer;

stateq = quantizer([q.StateWordLength q.StateFracLength]);
% The values are distributed over the range of the state format
zi = randquant(stateq, size(Hd.HiddenStates)); 

num = Hd.privNum;
den = Hd.privDen;
sv = Hd.privScaleValues;
issvnoteq2one = Hd.issvnoteq2one;
nsecs = Hd.nsections;

if length(issvnoteq2one) > nsecs + 1
    warning('dfilt:abstractsos:filter:ExtraScaleValues', ...
        'Cannot use more Scale Values than the number of sections plus one.');
    issvnoteq2one = issvnoteq2one(1:nsecs+1);
end

% Reach steady state
N = impzlength(Hd);
[y,zii] = df2sosfilter(q,num,den,sv,issvnoteq2one,x(1:N,:),zi);

% Look for limitcycles in steady state only
[y,zf,overflows] = df2sosfilter(q,num,den,sv,issvnoteq2one,x(N+1:end),zii,'limitcycle');

Hd.HiddenStates = zf;

% Return visible initial conditions only
[M,N] = size(zi);
zi = reshape(zi,2,M*N/2);

% [EOF]
