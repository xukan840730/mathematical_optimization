function ha = analogresp(h)
%ANALOGRESP   Compute analog response object.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:14:33 $

% Compute analog frequency
wp = cot(pi*h.Fpass/2);
ws = cot(pi*h.Fstop/2);

% Construct analog specs object
ha = fspecs.alppassfstop(h.FilterOrder,wp,ws,h.Apass);


% [EOF]
