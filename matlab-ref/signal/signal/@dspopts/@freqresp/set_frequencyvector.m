function fv = set_frequencyvector(this, fv)
%SET_FREQUENCYVECTOR   PreSet function for the 'frequencyvector' property.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2004/12/26 22:11:00 $

set(this, 'FrequencySpecification', 'FrequencyVector');

% if strcmpi(this.FrequencySpecification, 'NFFT')
%     siguddutils('readonlyerror', 'FrequencyVector', ...
%         'FrequencySpecification', 'FrequencyVector');
% end

set(this, 'privFrequencyVector', fv);

fv = [];

% [EOF]
