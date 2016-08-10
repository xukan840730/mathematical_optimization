function b = thisgenmcode(d)
%GENMCODE Generate M-Code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.2.6.2 $  $Date: 2004/04/13 00:09:10 $

h = d.ResponseTypeSpecs;
[params, values, descs, str, args] = genmcode(h, d);

[P,DENS] = getNumericSpecs(d);
params = {'Nb', 'Na', params{:}, 'P', 'dens'};
values = {getmcode(d, 'numOrder'), getmcode(d, 'denOrder'), values{:}, ...
        getmcode(d, P), getmcode(d, DENS)};
descs  = {'', '', descs{:}, 'P''th norm', ''};

IN = get(d,'initNum');
ID = get(d,'initDen');
if isempty(IN),
    in      = '';
	optargs = '';
else
    optargs = ', IN, ID';
    params  = {params{:}, 'IN', 'ID'};
    values  = {values{:}, sprintf('%d', IN), sprintf('%d', ID)};
    descs   = {descs{:}, '', ''};
end

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams(params, values, descs));
b.addcr(str, designdesc(d));
b.addcr('[b,a,err,sos_var,g] = iirlpnorm(%s);', ...
    sprintf('Nb, Na, %s, P, {dens}%s', args, optargs));
b.add('Hd                  = dfilt.df2sos(sos_var, g);');

% [EOF]
