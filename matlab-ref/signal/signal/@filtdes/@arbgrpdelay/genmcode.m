function b = genmcode(h, d)
%GENMCODE Generate M code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:01:36 $

params = {'N', 'F', 'E', 'G', 'W', 'R', 'P', 'dens'};
values = {getmcode(d, 'Order'), getmcode(d,'FrequencyVector'), ...
    getmcode(d,'FrequencyEdges'), getmcode(d,'GroupDelayVector'), ...
    getmcode(d,'WeightVector'), getmcode(d,'maxRadius'), ...
    getmcode(d, [get(d, 'initPnorm') get(d,'Pnorm')]), ...
    getmcode(d,'DensityFactor')};
descs  = {'', '', '', '', '', 'Maximum Pole Radius', 'P''th Norm', ''};

if isempty(d.InitDen),
    opt    = '';
    optstr = '';
else
    opt    = ', ID';
    params = {params{:}, 'ID'};
    values = {values{:}, getmcode(d, 'InitDen')};
    descs  = {descs{:}, 'Initial Denominator'};
end

fsstr = getfsstr(d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams(params, values, descs));
b.cr;
b.addcr('% Use the SOS matrix to avoid numerical round-off errors.');
b.addcr('[b, a, err, s] = iirgrpdelay(N, F%s, E%s, G, W, R, P, {dens}%s);', fsstr, fsstr, opt);
b.add('Hd             = dfilt.df2sos(s);');

% [EOF]
