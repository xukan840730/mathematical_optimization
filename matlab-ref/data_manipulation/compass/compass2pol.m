function theta_pol = compass2pol(theta)
%COMPASS2POL convert direction data (degN) into
% polar coordinates.
%   COMPASS2POL(THETA) convert the direction theta
%      (degree North) into polar coordinate.
%      note: theta is in degrees and between 0 and 360.
%
%   See also POL2CART
%

% Author: Arnaud Laurent
% Creation : March 20th 2009
% MATLAB version: R2007b
%

idx = find(theta>=0&theta<90);
theta_pol(idx,1) = abs(theta(idx) - 90);

idx = find(theta>=90&theta<=360);
theta_pol(idx,1) = abs(450 - theta(idx));

