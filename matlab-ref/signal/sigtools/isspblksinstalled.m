function [b, errstr, errid] = isspblksinstalled
%ISSPBLKSINSTALLED   Returns true if Signal Processing Blockset is installed.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2004/07/14 06:47:56 $

b = license('test', 'Signal_Blocks') && ~isempty(ver('dspblks'));

if b
    errstr = '';
    errid  = '';
else
    errstr = sprintf('%s\n%s', 'Signal Processing Blockset is not available.', ...
        'Make sure that it is installed and that a license is available.');
    errid  = 'noSPBlks';
end

% [EOF]
