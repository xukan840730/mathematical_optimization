function [y, t] = stepz(Hd, varargin)
%IMPZ Returns the impulse response

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.3.4.2 $  $Date: 2004/04/12 23:56:27 $

[y, t] = timeresp(Hd, @lclstepz, varargin{:});


% -------------------------------------------
function [y, t] = lclstepz(Hd, N, Fs)

if isa(Hd, 'qfilt'),
    
    if isempty(Fs), Fs = 1; end
    
    N = floor(N);
    
    t = (0:N-1)/Fs;
    
    % Multiply by .25 to avoid overflow
    x = .25*ones(N,1);
    y = filter(Hd,x)*4;
    
    varargout = {y, t};
else
    
    [y, t] = stepz(Hd, N, Fs);
    varargout = {y,t};
end

% [EOF]
