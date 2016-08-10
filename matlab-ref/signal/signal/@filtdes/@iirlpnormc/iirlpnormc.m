function d = iirlpnormc
%IIRLPNORMC  Constructor for this design method object.
%
%   Outputs:
%       d - Handle to the design method object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.3 $  $Date: 2002/04/15 00:43:16 $



d = filtdes.iirlpnormc;

% Call the super constructor
iirlpnorm_construct(d);

% Overwrite the tag
set(d,'Tag','IIR Constrained least P-th norm');








