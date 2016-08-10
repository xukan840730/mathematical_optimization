function notification_listener(hObj, eventData)
%NOTIFICATION_LISTENER Listener to the Notification event

%   Author(s): J. Schickler
%   Copyright 1988-2006 The MathWorks, Inc.
%   $Revision: 1.2.4.2 $  $Date: 2006/06/27 23:41:10 $

% Get the Notification type and all it's possible settings
NTypes = set(eventData, 'NotificationType');
NType  = get(eventData, 'NotificationType');

% Switch on the Notification type
switch NType
case NTypes{1}, % 'ErrorOccurred'
    error(hObj, 'WinTool Error', eventData.Data.ErrorString);
case NTypes{2}, % 'WarningOccurred'
    warning(hObj, 'WinTool Warning', eventData.Data.WarningString, ...
        eventData.Data.WarningID);
case NTypes{3}, % 'StatusChanged'
    % NO OP We are ignoring statuses for now.  See G 121740
case NTypes{4}, % 'FileDirty'
    % NO OP WINTool does not have files
otherwise
    error(hObj, ...
        ['Unhandled notification ''' NType ''' sent from ' class(eventData.Source) '.']);
end

% [EOF]
