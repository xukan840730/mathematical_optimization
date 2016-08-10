function varargout = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.4.4.6 $  $Date: 2005/12/22 19:02:08 $

% Set up design params
N = get(d,'order');

[Fpass, Apass] = getdesignspecs(h, d);

if nargout == 1,
    hfdesign = fdesign.lowpass('N,Fp,Ap', d.Order, Fpass, Apass);
    Hd       = cheby1(hfdesign);
    
    varargout = {Hd};
else

    [z,p,k] = cheby1(N,Apass,Fpass);
    varargout = {z,p,k};
end

% [EOF]
