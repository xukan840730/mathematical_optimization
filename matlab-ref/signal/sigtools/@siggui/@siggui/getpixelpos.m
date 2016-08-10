function pos = getpixelpos(this, field, varargin)
%GETPIXELPOS Get the position in pixel units.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/12/14 15:19:39 $

error(nargchk(2,inf,nargin,'struct'));

if ischar(field),
    field = this.Handles.(field);
    for indx = 1:length(varargin)
        if ischar(varargin{indx});
            field = field.(varargin{indx});
        else
            field = field(varargin{indx});
        end
    end
end

origUnits = get(field, 'Units');
set(field, 'Units', 'Pixels');
pos = get(field, 'Position');
if ~iscell(origUnits), origUnits = {origUnits}; end
for indx = 1:length(origUnits),
    set(field(indx), 'Units', origUnits{indx});
end

% [EOF]
