function nsecs = thisnsections(Hd)
%THISNSECTIONS Number of sections in a discrete filter.
%   THISNSECTIONS(Hd) returns the number of sections in a discrete
%   filter.
%
%   Example:
%       [b,a] = butter(7,.5);
%       Hd = sos(dfilt.df2t(b,a));
%       nsections(Hd)
% 
%   See also DFILT.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.8.3 $  $Date: 2005/06/16 08:18:04 $

% This should be private

% NSECTIONS is always 1 unless the class extends dfilt.multisection
nsecs = 1;

% [EOF]
