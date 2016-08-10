function Hd = design(h, d)
%DESIGN Design the filter

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:06:37 $

N = get(d, 'Order');

Fc = get(d, 'Fc');

Upper = [get(d, 'DstopUpper') 1+get(d, 'DpassUpper')];
Lower = [-get(d, 'DstopLower') 1-get(d, 'DpassLower')];

b = fircls(N, [0 Fc 1], [0 1], Upper, Lower);

Hd = dfilt.dffir(b);

% [EOF]
