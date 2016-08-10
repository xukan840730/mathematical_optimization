function centerdc = set_centerdc(this, centerdc)
%SET_CENTERDC   PreSet function for the 'CenterDC' property.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2004/01/25 23:06:59 $

this.privcenterdc = centerdc;

% Force to use the entire nyquist interval when appropriate
if centerdc,
    fullnyq(this);
end

% Don't duplicate
centerdc = [];


% [EOF]
