function b = genmcode(h, d)
%GENMCODE Generate M-Code

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.
%   $Revision: 1.1.4.6.32.1 $  $Date: 2009/12/17 13:59:39 $

b = sigcodegen.mcodebuffer;

[p,v] = abstract_genmcode(h, d);
p{end+1} = 'match';
v{end+1} = sprintf('''%s''', get(d, 'MatchExactly'));

b.addcr(b.formatparams(p,v));
b.cr;
b.addcr(designdesc(d));
b.addcr('h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop%s);', getfsinput(d));
b.add('Hd = design(h, ''cheby1'', ''MatchExactly'', match);');

% [EOF]