function Hd = design(h, d)
%DESIGN Design the filter

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:06:42 $

N = get(d, 'Order');

Fc = get(d, 'Fc');

Upper = [1+get(d, 'DpassUpper') get(d, 'DstopUpper')];
Lower = [1-get(d, 'DpassLower') -get(d, 'DstopLower')];

b = fircls(N, [0 Fc 1], [1 0], Upper, Lower);

Hd = dfilt.dffir(b);

% [EOF]
