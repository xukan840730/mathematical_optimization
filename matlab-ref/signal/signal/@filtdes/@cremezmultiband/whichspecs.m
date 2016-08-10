function ws = whichspecs(h)
%WHICHSPECS

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.1.8.1 $  $Date: 2003/03/02 10:15:44 $

ws = ft_whichspecs(h);

ws(1).defval = [-24000 -12000 -9600 16800 19200 24000];
ws(2).defval = [0 0 1 2 0 0];
ws(3).defval = [1 1 1];

% [EOF]
