function normalizefreq(this,normFlag, Fs)
%NORMALIZEFREQ   Normalize/un-normalize the frequency of the data object.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/04/09 19:05:14 $

% This method doesn't do anything too interesting, but its overloaded by
% its subclasses.
if nargin > 1
    this.privNormalizedFrequency = normFlag;
    if nargin > 2
        this.Fs = Fs;
    end
end

% [EOF]
