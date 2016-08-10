function setsosprivvals(this, privvals)
%SETSOSPRIVVALS   

%   Author(s): R. Losada
%   Copyright 2003 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2003/12/06 16:00:45 $

% Get names of private properties we want to copy
pnames = sosprivnames(this);

set(this,pnames,privvals.this);

% Set parent properties last. We need nsections to be correct
setdfiltprivvals(this, privvals.parent);



% [EOF]
