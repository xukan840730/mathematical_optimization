function examples = getexamples(this)
%GETEXAMPLES   Get the examples.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/10/23 18:50:20 $

examples = {{ ...
    'Compare passband and stopband MatchExactly.', ...
    'h     = fdesign.bandpass(''Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2'');', ...
    'Hd    = design(h, ''butter'', ''MatchExactly'', ''passband'');', ...
    'Hd(2) = design(h, ''butter'', ''MatchExactly'', ''stopband'');', ...
    '', ...
    '% Compare the passband edges in FVTool.', ...
    'fvtool(Hd);', ...
    'axis([.44 .56 -2 0]);'}};

% [EOF]