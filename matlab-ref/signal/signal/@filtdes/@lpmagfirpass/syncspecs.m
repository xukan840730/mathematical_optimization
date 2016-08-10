function props = syncspecs(h,d)
%SYNCSPECS Properties to sync.
%
%   Inputs:
%       h - handle to this object
%		d - handle to container objects
%
%   Outputs:
%       props - cell array, with list of properties to sync from
%               and list of properties to sync to.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:09:55 $

magUnits = get(d,'magUnits');
magUnitsOpts = set(d,'magUnits');

switch magUnits,
case magUnitsOpts{1}, % 'dB'
	props = {'Apass'};
case magUnitsOpts{2}, % 'Linear'
	props = {'Dpass'};
end	

% [EOF]
