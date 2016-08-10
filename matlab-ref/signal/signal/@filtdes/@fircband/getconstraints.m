function varargout = getconstraints(hObj)
%CONVERTCEM Convert the Constrained Error Magnitudes

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:04:34 $

cb = get(hObj, 'ConstrainedBands');

if isempty(cb),
    cem = [];
else
    n = nbands(hObj.ResponseTypeSpecs, hObj);
    if any(cb > n),
        error(generatemsgid('BandOutOfBounds'), ...
            'Cannot constrain the %s band.  There are only %d bands.', ...
            getplace(max(cb)), n);
    end
    if length(cb) >= n,
        error(generatemsgid('TooManyConstraints'), ...
            'FIRCBAND cannot constrain all bands.');
    end
    cem = repmat({'w'}, 1, n);
    for indx = 1:length(cb)
        cem{cb(indx)} = 'c';
    end
end

if nargout == 1,
    varargout = {cem};
else
    if ~isempty(cem), cem = genmcodeutils('formatcellstr', cem); end
    varargout = {'CEMs', cem, 'Constrained Error Magnitudes'};
end

% --------------------
function p = getplace(n)

if n > 3 && n < 21
    p = sprintf('%dth', n);
else
    switch rem(n,10)
        case 1
            p = sprintf('%dst', n);
        case 2
            p = sprintf('%dnd', n);
        case 3
            p = sprintf('%drd', n);
        otherwise
            p = sprintf('%dth', n);
    end
end

% [EOF]
