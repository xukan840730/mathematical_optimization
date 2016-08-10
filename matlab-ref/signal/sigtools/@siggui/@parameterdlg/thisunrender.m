function thisunrender(this)
%THISUNRENDER Unrender for the parameter dialog

%   Author(s): J. Schickler
%   Copyright 1988-2008 The MathWorks, Inc.
%   $Revision: 1.2.4.5 $  $Date: 2009/01/05 18:00:54 $

delete(this.BaseListeners);

hFig = get(this, 'FigureHandle');
if ~isempty(hFig) && ishghandle(hFig),
    delete(hFig);
end

% Not sure why this was here but it was causing problems.
% % Reset the parameters.
% hPrm = get(this, 'Parameters');
% for indx = 1:length(hPrm)
%     send(hPrm(1), 'UserModified', sigdatatypes.sigeventdata(hPrm(1), ...
%     'UserModified', get(hPrm(1), 'Value')));
% end

% [EOF]
