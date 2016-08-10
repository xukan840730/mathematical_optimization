function b = genmcode(h, d)
%GENMCODE Generate M-Code

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.
%   $Revision: 1.1.4.6.32.1 $  $Date: 2009/12/17 13:59:17 $

b = sigcodegen.mcodebuffer;

[p, v] = abstract_genmcode(h, d);
p{end+1} = 'match';
v{end+1} = sprintf('''%s''', get(d, 'MatchExactly'));

b.addcr(b.formatparams(p, v));
b.cr;
b.addcr(designdesc(d));
b.addcr('h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, Apass2%s);', getfsinput(d));
b.add('Hd = design(h, ''cheby2'', ''MatchExactly'', match);');

% [EOF]
