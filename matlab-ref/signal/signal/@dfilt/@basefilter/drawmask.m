function varargout = drawmask(this, varargin)
%DRAWMASK   Draw the mask for the filter specifications.

%   Author(s): J. Schickler
%   Copyright 2005-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2007/12/14 15:07:25 $

hd = getfdesign(this);

if isempty(hd)
    warning(generatemsgid('Empty'),'There are no specifications to draw');
    return;
end

% Make sure that we scale the mask if the filter has any interpolation.
if nargin > 1 && strcmpi(varargin{end}, 'normalize')
    normalize     = true;
    scalefactor   = 1;
    varargin(end) = [];
else
    normalize   = false;
    scalefactor = nominalgain(this);
end

hm = getfmethod(this);

[x, y, h, fcns] = drawmask(hd, hm, varargin{:});

if scalefactor == 1
    if normalize
        ydata = get(h, 'YData');
        maxy  = max(ydata);

        if strcmpi(lower(fcns.getunits()), 'db')
            if ~isnan(maxy) && isfinite(maxy)
                ydata = ydata-maxy;
            end
        else
            if maxy ~= 0 && ~isnan(maxy) && isfinite(maxy)
                ydata = ydata./maxy;
            end
        end
        set(h, 'YData', ydata);
    end
else

    ydata = get(h, 'YData');
    switch lower(fcns.getunits())
        case 'db'
            ydata = ydata+db(scalefactor);
        case 'squared'
            ydata = ydata*scalefactor^2;
        otherwise
            ydata = ydata*scalefactor;
    end

    set(h, 'YData', ydata);
end

if nargout == 2
    varargout = {x,y};
else
    varargout = {h};
end


% [EOF]
