function this = delay(lat)
%DELAY Integer delay.
%   Hd = DFILT.DELAY(D) constructs a discrete-time integer delay object
%   with a latency of D.
%
%   % EXAMPLE
%   Hd = dfilt.delay(2)
%
%   See also DFILT/STRUCTURES   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2005/12/22 18:57:16 $

this = dfilt.delay;
this.filterquantizer = dfilt.filterquantizer;
this.FilterStructure = 'Delay';

if nargin>=1
  this.Latency = lat;
end

% [EOF]
