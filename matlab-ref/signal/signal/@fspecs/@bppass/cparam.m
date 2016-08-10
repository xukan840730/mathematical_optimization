function c = cparam(h)
%CPARAM   

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:12:04 $

F1 = h.Fpass1;
F2 = h.Fpass2;
c = computecparam(h,F1,F2);


% [EOF]
