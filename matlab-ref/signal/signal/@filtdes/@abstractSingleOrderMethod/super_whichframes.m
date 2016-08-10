function fr = super_whichframes(h)
%WHICHFRAMES  Return constructors of frames needed for FDATool.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.5 $  $Date: 2002/04/15 00:29:14 $

% Call the base method
fr = dm_whichframes(h);

% Get frames needed by filter order
fr(end+1).constructor = 'siggui.filterorder';
fr(end).setops      = {'isMinOrd',0};


% Show no options by default
fr(end+1).constructor = 'siggui.textOptionsFrame';
fr(end).setops      = {};
