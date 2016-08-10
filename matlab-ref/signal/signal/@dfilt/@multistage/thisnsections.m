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
%   $Revision: 1.1.6.1 $  $Date: 2004/12/26 22:08:59 $

% This should be private

nsecs = sum(nsections(Hd.Stage));

% [EOF]
