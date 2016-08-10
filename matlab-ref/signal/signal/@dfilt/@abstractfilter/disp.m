function disp(this)
%DISP Object display.
  
%   Author: V. Pellissier
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.4.4.4 $  $Date: 2005/06/30 17:33:09 $

if length(this) > 1
    vectordisp(this);
    return;
end

fn = fieldnames(this);
N = length(fn);
% Reorder the fields. NumSamplesProcessed, ResetStates and States in
% the end.
fn = fn([3, 5:N, 1, 4, 2]);

siguddutils('dispstr', this, fn);

% [EOF]
