%OPTSIMINIT Sets up necessary data files for optimization of optsim.
% Define Optimization parameters

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/08/03 21:32:01 $


% Define Tunable Variable initial values
Kp = 0.63;
Ki = 0.0504;
Kd = 1.9688;

% Define Uncertain Variable initial values
a2 = 43;
a1 = 3;

disp('Done initializing optsim.')
% end optsiminit
