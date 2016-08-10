function theta_comp = pol2compass(theta)
%POL2COMPASS convert polar direction data (deg) into
% degN reference.
%   POL2COMPASS(THETA) convert the direction theta
%      (polar) into degree North coordinate.
%
%   See also POL2CART
%

% Author: Arnaud Laurent
% Creation : March 20th 2009
% MATLAB version: R2007b
%

idx = find(theta<0);
theta(idx) = 360 + theta(idx);

idx = find(theta>=0&theta<90);
theta_comp(idx,1) = abs(theta(idx) - 90);

idx = find(theta>=90&theta<=360);
theta_comp(idx,1) = abs(450 - theta(idx));
