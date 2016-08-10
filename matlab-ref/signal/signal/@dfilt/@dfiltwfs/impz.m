function [y, t] = impz(Hd, varargin)
%IMPZ Returns the impulse response

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.2.4.3 $  $Date: 2009/07/27 20:29:09 $

[y, t] = timeresp(Hd, @lclimpz, varargin{:});

% -----------------------------------------------------------
function [y, t] = lclimpz(G, N, Fs)

[y, t] = impz(G, N, Fs);

% [EOF]
