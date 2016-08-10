function [params, values, descs, str, args] = genmcode(h, d)
%GENMCODE Generate M code

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2002/10/04 18:13:52 $

fsstr = getfsstr(d);

params = {'F', 'E', 'A', 'W'};
values = {getmcode(d,'FrequencyVector'), getmcode(d,'FrequencyEdges'), ...
        getmcode(d,'MagnitudeVector'), getmcode(d,'WeightVector')};
descs  = {'', '', '', ''};

str = '';

args = sprintf('F%s, E%s, A, W', fsstr, fsstr);

% [EOF]
