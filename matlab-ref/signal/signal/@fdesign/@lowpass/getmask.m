function [F, A] = getmask(this, fcns, rcf, specs)
%GETMASK   Get the mask.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2007/06/14 05:23:16 $

% If the specs were not passed in or are [], use the design specifications.
if nargin < 4 || isempty(specs)
    specs = getspecs(this.CurrentSpecs);
    fpass = specs.Fpass;
    fstop = specs.Fstop;
else
    fpass = specs.Fpass;
    fstop = specs.Fstop;
    if ~specs.NormalizedFrequency
        fpass = fpass/specs.Fs*2;
        fstop = fstop/specs.Fs*2;
    end
end

if isempty(fpass)
    fpass = NaN;
end
if isempty(fstop)
    fstop = NaN;
end

% The frequency vector is always the same.
F = [1 fpass fpass 0 0 fstop fstop 1]*fcns.getfs()/2;

A = fcns.gethighlow(specs);

% [EOF]
