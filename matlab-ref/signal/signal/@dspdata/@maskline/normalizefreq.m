function normalizefreq(this, newvalue, fs)
%NORMALIZEFREQ   Normalize the frequency and the vector.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/12/14 15:10:23 $

error(nargchk(3,3,nargin,'struct'));

% If the new value == the old value return early.
if newvalue == this.NormalizedFrequency
    return;
end

if newvalue
    % Going to normalized.
    this.FrequencyVector = this.FrequencyVector./(fs/2);
else
    % Turning off normalized.
    this.FrequencyVector = this.FrequencyVector.*(fs/2);
end

this.NormalizedFrequency = newvalue;

% [EOF]
