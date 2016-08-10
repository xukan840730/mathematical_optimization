function p = propstoadd(this)
%PROPSTOADD   Return the properties to add to the parent object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/01/20 15:35:43 $

p = fieldnames(this);
p = {p{1:4},p{8:9},p{6},p{5},p{7},p{10:end}};

p(1) = []; % All but the responsetype.

% [EOF]
