function varargout = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.4.4.5 $  $Date: 2005/12/22 19:02:20 $

[Fpass, Fstop, Apass, Astop] = getdesignspecs(h,d);

if nargout == 1,
    hfdesign = fdesign.lowpass(Fpass, Fstop, Apass, Astop);
    Hd       = ellip(hfdesign, 'MatchExactly', d.MatchExactly);
    
    varargout = {Hd};
else
    N = ellipord(Fpass,Fstop,Apass,Astop);

    [z,p,k] = ellip(N,Apass,Astop,Fpass);
    varargout = {z,p,k};
end

% [EOF]
