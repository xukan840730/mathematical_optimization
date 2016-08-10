function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.2.4.1 $  $Date: 2003/03/02 10:23:20 $

args = getarguments(h, d);
dens = get(d,'DensityFactor');
N    = get(d, 'Order');

% Design a type 4 if order is odd
if rem(N, 2),
    flag = {'Hilbert'};
else
    flag = {};
end

b = remez(N, args{:}, {dens}, flag{:});

% Construct object
Hd = dfilt.dffir(b);

% [EOF]
