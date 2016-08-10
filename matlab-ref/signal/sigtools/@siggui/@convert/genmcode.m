function str = genmcode(hObj)
%GENMCODE Generate M-code

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2002/10/15 18:07:03 $

% Get the new filter structure name for the filter object.
structure = getconstructor(hObj);

str = sprintf('%s\n%s', ...
    sprintf('%% Convert the filter to the %s structure.', get(hObj, 'TargetStructure')), ...
    sprintf('Hd = convert(Hd, ''%s'');', structure));

% [EOF]
