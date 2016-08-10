function HdB = convert2db(this,H)
%CONVERT2DB   Convert input response to db values.

%   Author(s): P. Pacheco
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2004/12/26 22:10:04 $

ws = warning; % Cache warning state
warning off   % Avoid "Log of zero" warnings
HdB = db(H,'power');  % Call the Convert to decibels engine
warning(ws);  % Reset warning state

% [EOF]
