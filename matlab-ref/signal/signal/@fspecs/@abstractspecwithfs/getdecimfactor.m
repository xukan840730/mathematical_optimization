function decimfactor = getdecimfactor(this)
%GETDECIMFACTOR   Get the decimfactor.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2005/04/04 17:00:10 $

if isprop(this, 'DecimationFactor')
    decimfactor = get(this, 'DecimationFactor');
elseif isprop(this, 'Band')
    decimfactor = get(this, 'Band');
else
    % We should probably error here instead, but this will help the
    % HALFBAND case, which does not have a property for the band because it
    % is always 2.
    decimfactor = 2;
end

% [EOF]
