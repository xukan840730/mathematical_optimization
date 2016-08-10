function info = iir_qtoolinfo(this)
%IIR_QTOOLINFO   Return the IIR specific info.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2004/01/25 23:04:37 $

info.normalize = 'numerator';

info.coeff.setops  = {'Name', 'Coefficient', ...
    'FracLabels', {'Numerator', 'Denominator'}};
info.coeff.syncops = {'Num', 'Den'};

info.product.setops  = {'FracLabels', {'Num.', 'Den.'}};
info.product.syncops = {'NumProd', 'DenProd'};

info.accum.setops  = {'FracLabels', {'Num.', 'Den.'}};
info.accum.syncops = {'NumAccum', 'DenAccum'};

% [EOF]
