function orderMode_update(h,ft)
%ORDERMODE_UPDATE Update the order mode property.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.3.4.1 $  $Date: 2008/05/31 23:26:47 $

constr = findConstr(h,ft,'minimum');

if isempty(constr),
    % Disable orderMode
    enabdynprop(h,'orderMode','off');
else
    % Get window
    win = get(h,'Window');
    
    % If window is not kaiser, disable ordermode 
    if ~strcmpi(win,'kaiser'),
        enabdynprop(h,'orderMode','off');
    else
        % Enable orderMode
        enabdynprop(h,'orderMode','on');      
    end
end

% Update the window parameter
update_winparam(h);

