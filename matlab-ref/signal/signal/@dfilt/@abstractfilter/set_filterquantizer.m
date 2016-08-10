function newfq = set_filterquantizer(this, newfq)
%SET_FILTERQUANTIZER   PreSet function for the 'filterquantizer' property.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2007/12/14 15:06:59 $

oldfq = get(this, 'filterquantizer');

set(this, 'privfilterquantizer', newfq);

% Remove current dynamic properties
rmprops(this, oldfq);

set_ncoeffs(newfq, naddp1(this));

try
    % Quantize the coefficients
    quantizecoeffs(this);
catch ME
    
    % Get the old quantizer because the new one errors.
    set(this, 'privfilterquantizer', oldfq);

    newfq = oldfq;
    quantizecoeffs(this);
    
    addprops(this, newfq);
    throwAsCaller(ME);
end

addprops(this, newfq);

validatestates(this);

% Store nothing, its stored in private.
newfq = [];

% [EOF]
