function h = reffilter(this)
%DOUBLE   Returns the double representation of the filter object.
%   See help in dfilt/reffilter.

%   Author(s): R. Losada
%   Copyright 2003-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/04/21 04:35:15 $

% Get array of contained dfilts
secs = this.Stage;

for n = 1:length(secs),
    newsecs(n) = reffilter(secs(n));
end

h = feval(str2func(class(this)),newsecs(:));

% Set fdesign/fmethod objs in new obj
setfdesign(h, getfdesign(this));
setfmethod(h, getfmethod(this));



% [EOF]
