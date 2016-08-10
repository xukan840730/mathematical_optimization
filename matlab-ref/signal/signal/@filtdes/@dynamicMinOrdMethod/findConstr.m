function consStr = findConstr(h,ft,orderMode)
%FINDCONSTR Find the appropriate constructor.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 2002/04/15 00:36:07 $

s = get(h,'availableTypes');

indx = findConstrIndx(h,ft);

if nargin < 3,
    orderMode = 'specify';
end

% Return constructor
consStr = getfield(s(indx).construct,orderMode);