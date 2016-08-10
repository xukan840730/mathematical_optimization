function examples = getexamples(this)
%GETEXAMPLES   Get the examples.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2005/06/30 17:37:59 $

examples = {{ ...
    'Design a bandpass Equiripple filter with a weighted passband.', ...
    'h  = fdesign.bandpass(''N,Fst1,Fp1,Fp2,Fst2'', 30);', ...
    'Hd = design(h, ''equiripple'');', ...
    '', ...
    '% Weight the passband more than the stopbands.', ...
    'Hd(2) = design(h, ''equiripple'', ''Wpass'', 20);'....
    '', ...
    '% Compare the two designs.', ...
    'fvtool(Hd)'}};

% [EOF]
