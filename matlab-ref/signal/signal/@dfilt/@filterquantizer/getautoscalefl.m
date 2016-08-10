function fl = getautoscalefl(this,s,signed,wl)
%GETAUTOSCALEFL   Get the autoscalefl.

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2006/06/27 23:34:31 $

if isfinite(s.Min) && s.Min~=realmax,
    f1 = fi(s.Min,signed,wl);
else
    f1.FractionLength = wl-1;
end
if isfinite(s.Max) && s.Max~=-realmax,
    f2 = fi(s.Max,signed,wl);
else
    f2.FractionLength = wl-1;
end
fl = min(f1.FractionLength,f2.FractionLength);

% [EOF]
