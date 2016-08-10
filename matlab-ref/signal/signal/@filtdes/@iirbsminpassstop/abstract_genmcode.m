function varargout = abstract_genmcode(h, d)
%ABSTRACT_GENMCODE Returns the F and A M-code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2004/04/13 00:09:05 $

[Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, Apass2] = getdesignspecs(h,d);

p = {'Fpass1', 'Fstop1', 'Fstop2', 'Fpass2', 'Apass1', 'Astop', 'Apass2'};
v = {getmcode(d, Fpass1), getmcode(d, Fstop1), getmcode(d, Fstop2), ...
        getmcode(d, Fpass2), getmcode(d, Apass1), getmcode(d, Astop), getmcode(d, Apass2)};

if nargout == 1
    b = sigcodegen.mcodebuffer;

    b.addcr(b.formatparams(p, v));
    varargout = {b};
else
    varargout = {p, v};
end

% [EOF]
