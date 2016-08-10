function varargout = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.4.4.8 $  $Date: 2007/12/14 15:12:23 $

[Fpass1, Fpass2, Apass] = getdesignspecs(h, d);

if nargout == 1
    hfdesign = fdesign.bandstop('N,Fp1,Fp2,Ap', d.Order, Fpass1, Fpass2, Apass);
    Hd       = cheby1(hfdesign);
    
    varargout = {Hd};
else
    % Set up design params
    N = get(d,'order');

    if rem(N,2),
        error(generatemsgid('MustBeEven'),'Bandstop designs must be of even order.');
    end

    F = [Fpass1 Fpass2];

    [z,p,k] = cheby1(N/2,Apass,F,'stop');

    varargout = {z,p,k};
end

% [EOF]
