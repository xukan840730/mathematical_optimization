function dBcase(h,d)
%DBCASE Handle the dB case.
%
% This should be a private method.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:02:57 $

convertmag(h,d,...
    {'Dpass1', 'Dpass2'},...
    {'Apass1', 'Apass2'},...
    {'pass', 'pass'},...
    'todb');

% [EOF]
