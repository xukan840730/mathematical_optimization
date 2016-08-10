function [F, A, W, args] = getarguments(h, d)
%GETARGUMENTS Return the design method arguments

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:04:51 $

F = [0 get(d, 'Fstop1') get(d, 'Fpass1') get(d, 'Fpass2') get(d, 'Fstop2') 1];
A = [0 0 1 1 0 0];

mu = get(d, 'magUnits'); set(d, 'magUnits', 'linear');
W  = [1 get(d, 'Dpass') get(d, 'Dstop2')]; set(d, 'magUnits', mu);

args = {};

% [EOF]
