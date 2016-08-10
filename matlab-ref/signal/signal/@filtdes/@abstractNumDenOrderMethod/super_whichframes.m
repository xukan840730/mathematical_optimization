function fr = super_whichframes(h)
%WHICHFRAMES  Return constructors of frames needed for FDATool.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 2002/03/28 17:10:40 $

% Call the base method
fr = dm_whichframes(h);

% Get frames needed by filter order
fr(end+1) = whichframes(get(h,'numDenFilterOrderObj'));
fr(end).setops      = {};


% Show no options by default
fr(end+1).constructor = 'siggui.textOptionsFrame';
fr(end).setops      = {};
