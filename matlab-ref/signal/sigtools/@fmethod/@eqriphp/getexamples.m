function examples = getexamples(this)
%GETEXAMPLES   Get the examples.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2005/06/30 17:38:16 $

examples = {{ ...
    'Design a highpass FIR Least-squares filter with a weighted stopband.', ...
    'h  = fdesign.highpass(''N,Fst,Fp'', 30);', ...
    'Hd = design(h, ''equiripple'');', ...
    '', ...
    '% Weight the stopband more than the passbands.', ...
    'Hd(2) = design(h, ''equiripple'', ''Wstop'', 20);'....
    '', ...
    '% Compare the two designs.', ...
    'fvtool(Hd)'}};

% [EOF]
