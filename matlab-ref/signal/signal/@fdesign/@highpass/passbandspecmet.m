function flag = passbandspecmet(Hf,Hd,ng)
%PASSBANDSPECMET Check whether passband response is within spec.
%   This should be a private method.

%   Author(s): R. Losada
%   Copyright 2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/04/21 04:35:52 $

% Gather specs for measuring
mi = measureinfo(Hf);

if ~isempty(mi.Apass) && ~isempty(mi.Fpass),
    % Note: Apass could be empty even with Fpass defined (e.g. constrained
    % equiripple)
    w = linspace(pi*mi.Fpass,pi,512);
    H = freqz(Hd,w);
    absH = abs(H);
    
    % Convert passband peak-to-peak ripple to linear units
    Rpass = convertmagunits(mi.Apass,'db','linear','pass');
    
    flag = ~any((absH > ng+Rpass) | (absH < ng-Rpass));
else
    flag = true;
end



% [EOF]
