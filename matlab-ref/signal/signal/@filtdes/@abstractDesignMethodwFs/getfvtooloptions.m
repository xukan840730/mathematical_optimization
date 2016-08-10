function opts = getfvtooloptions(d)
%GETFVTOOLOPTIONS Returns the options sent to FVTool from the design method

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.4 $  $Date: 2002/11/21 15:39:01 $

opts = {'Analysis', 'magnitude'};

analysis = 'magnitude';

if isprop(d, 'magUnits'),
    switch lower(get(d, 'magUnits'))
    case 'db'
        opts = {opts{:}, 'MagnitudeDisplay', 'magnitude (db)'};
    case 'linear'
        opts = {opts{:}, 'MagnitudeDisplay', 'magnitude'};
    case 'squared'
        opts = {opts{:}, 'MagnitudeDisplay', 'magnitude squared'};
    end
else
    opts = {opts{:}, 'MagnitudeDisplay', 'magnitude'};
end

% [EOF]
