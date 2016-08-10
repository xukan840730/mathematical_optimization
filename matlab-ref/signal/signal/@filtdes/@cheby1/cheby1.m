function d = cheby1
%CHEBY1  Constructor for this design method object.
%
%   Outputs:
%       d - Handle to the design method object

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.3.4.3 $  $Date: 2004/04/13 00:03:57 $


d = filtdes.cheby1;

% Call super's constructor
classiciir_construct(d);

% Set the tag
set(d,'Tag','Chebyshev Type I', 'MatchExactly', 'passband');

