function n = nmult(this,optimones,optimnegones)
%NMULT Returns the number of multipliers  

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2005/08/20 13:25:41 $

n = 0;
for k=1:length(this.Stage)
    n = n + nmult(this.Stage(k),optimones,optimnegones);
end



% [EOF]
