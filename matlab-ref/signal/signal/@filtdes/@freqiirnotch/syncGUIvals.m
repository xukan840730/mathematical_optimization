function syncGUIvals(h, d, arrayh)
%SYNCGUIVALS Sync the values from the GUI

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2003/03/02 10:18:32 $

% Get frame to sync from
frames = whichframes(h);

% Find the right frame
fr = find(arrayh,'-class',frames.constructor);

% Make sure to sync the transition mode first
set(d, 'TransitionMode', get(fr,'TransitionMode'));

% Sync bandwidth or rolloff
set(d, get(d,'TransitionMode'), evaluatevars(get(fr,get(fr,'TransitionMode'))));

% Allow the subclasses to sync too
thissyncGUIvals(h,d,fr);

% [EOF]
