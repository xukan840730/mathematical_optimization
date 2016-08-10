function wf = whichframes(this)
%WHICHFRAMES Returns the frames for FDATool to use.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2004/12/26 22:12:21 $

wf = dmom_whichframes(this);

indx = find(strcmpi('siggui.textOptionsFrame', {wf.constructor}));

wf(indx).constructor = 'siggui.gremezoptsframe';
wf(indx).setops      = {'DisabledProps', disabledprops(this)};

if ~any(strcmpi(this.ResponseType, {'bandstop', 'highpass'}))
    indx = find(strcmpi('siggui.filterorder', {wf.constructor}));
    
    % Replace the filterorder frame with the gremez specific frame.
    wf(indx).constructor = 'siggui.gremezfilterorder';
end

% [EOF]
