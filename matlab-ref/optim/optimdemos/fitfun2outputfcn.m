function stop = fitfun2outputfcn(x,optimvalues,state,Data,plothandle)
%FITFUN2OUTPUTFCN Output function for FITFUN2 in DATDEMO.

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/08/03 21:31:50 $

stop = false;
if strcmp(state,'iter')
  % Call f to get the current yEstimate since that is not available in
  % OPTIMVALUES.
  [fval, yEstimate] = fitfun2(x,Data); 
  set(plothandle, 'ydata', yEstimate);
  drawnow; % Draws current graph now
end

