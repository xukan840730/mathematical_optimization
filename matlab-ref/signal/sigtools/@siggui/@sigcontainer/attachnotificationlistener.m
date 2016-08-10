function attachnotificationlistener(hParent)
%ATTACHNOTIFICATIONLISTENER

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2002/06/18 13:20:38 $

hAllChildren = allchild(hParent);

% Add a listener to a local function.  Creating function handles for
% external mfiles is very slow.  Local functions is much faster.
hListener = handle.listener(hAllChildren, 'Notification', @lclnotification_listener);
set(hListener, 'CallbackTarget', hParent);

set(hParent, 'NotificationListener', hListener);

% -----------------------------------------------------------
function lclnotification_listener(hObj, eventData, varargin)

notification_listener(hObj, eventData, varargin{:});

% [EOF]
