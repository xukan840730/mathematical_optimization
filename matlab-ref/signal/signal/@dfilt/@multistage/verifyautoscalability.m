function verifyautoscalability(this)
%VERIFYAUTOSCALABILITY   

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/06/11 17:30:53 $

for k=1:length(this.Stage)
  verifyautoscalability(this.Stage(k));
end

% [EOF]
