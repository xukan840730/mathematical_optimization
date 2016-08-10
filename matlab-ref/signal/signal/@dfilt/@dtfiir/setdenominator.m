function den = setdenominator(Hd, den)
%SETDENOMINATOR Overloaded set on the Denominator property.
  
%   Author: V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.3.4.5 $  $Date: 2007/04/09 19:05:10 $

if den(1)==0,
    error(generatemsgid('LeadDenCoeff'), ...
        'The leading coefficient of the denominator can''t be equal to 0.');
end

den = set_coeffs(Hd, den);

ncoeffs = Hd.ncoeffs;
oldlength=0;
if length(ncoeffs)==2, oldlength = ncoeffs(2); end
newlength = length(den);
ncoeffs(2) = newlength;
Hd.ncoeffs = ncoeffs;

% Update the ncoeffs property of the plugin.
set_ncoeffs(Hd.filterquantizer, ncoeffs);

if newlength~=oldlength,
    reset(Hd);
end

% Store numerator as reference
set(Hd,'refden',den);

% Quantize the coefficients
quantizecoeffs(Hd);

% Hold an empty to not duplicate storage
den = [];
