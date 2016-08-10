function [fsstr, fs] = getfsstr(d)
%GETFSSTR Returns '/(Fs/2)' or ''

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2002/10/04 18:12:55 $

if isnormalized(d),
    fsstr = '';
    fs    = '1';
else
    fsstr = '/(Fs/2)';
    fs    = 'Fs/2';
end

% [EOF]
