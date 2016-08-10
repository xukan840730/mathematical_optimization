function mv = get_magnitudevector(this, mv)
%GET_MAGNITUDEVECTOR   PreGet function for the 'magnitudevector' property.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2004/12/26 22:10:21 $

mv = get(this, 'privMagnitudeVector');

switch lower(this.MagnitudeUnits)
    case 'db'
        mv = db(mv);
    case 'linear'
        % NO OP.
    case 'squared'
        mv = mv.^2;
end

% [EOF]
