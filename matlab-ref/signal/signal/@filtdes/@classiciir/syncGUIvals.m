function syncGUIvals(this, arrayh)
%SYNCGUIVALS   Sync the values from the GUI.

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.3.54.1 $  $Date: 2009/12/17 13:59:21 $

super_syncGUIvals(this, arrayh)

if ~isspecify(this),
    wf = whichframes(this);
    h  = find(arrayh, '-class', wf(end).constructor);

    % We only want to do this if filterdesign exists.
    if ~isempty(h)
        set(this, 'MatchExactly', get(h, 'MatchExactly'));
    end
end

% [EOF]
