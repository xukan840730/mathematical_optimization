function ha = analogresp(h)
%ANALOGRESP   Compute analog response object.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:12:10 $

% Compute analog frequency
c = cparam(h);
wp1 = (c-cos(pi*h.Fpass1))/sin(pi*h.Fpass1);
wp2 = (c-cos(pi*h.Fpass2))/sin(pi*h.Fpass2);
wp = min(abs([wp1,wp2]));

% Construct analog specs object
ha = fspecs.alppassastop(h.FilterOrder,wp,h.Apass,max(h.Astop1,h.Astop2));


% [EOF]
