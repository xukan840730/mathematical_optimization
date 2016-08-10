function iirlpnorm_construct(d)
%IIRLPNORM  Constructor for this design method object.
%

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.3 $  $Date: 2002/04/15 00:42:56 $

% Call super's constructor
lpnorm_construct(d);

% Construct a numDenFilterOrder object and store it
h = filtdes.numDenFilterOrder;
set(d,'numDenFilterOrderObj',h);

% Set the tag
set(d,'Tag','IIR Least P-th norm');

