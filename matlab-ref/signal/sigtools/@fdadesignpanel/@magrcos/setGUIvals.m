function setGUIvals(this, eventData) %#ok
%SETGUIVALS Set the values in the GUI

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.2.4.4.4.1 $  $Date: 2009/12/14 20:10:57 $

h = findhandle(this, whichframes(this));

if ~isempty(h),
    
    if strcmpi(this.DesignType, 'normal')
        str = {'The attenuation at cutoff', 'frequencies is fixed at 6 dB', ...
            '(half the passband gain)'};
    else
        str = {'The attenuation at cutoff', 'frequencies is fixed at 3 dB', ...
            '(half the passband power)'};
    end
    set(h, 'Comment', str, ...
    	'AllOptions', set(this, 'DesignType'), ...
        'currentSelection', get(this, 'DesignType'));
end

% [EOF]
