function specs = whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:03:02 $

% Call super's method
specs = mf_whichspecs(h);

% Prop name, data type, default value, listener callback
specs(end+1) = cell2struct({'Apass1','udouble',1,[],'magspec'},specfields(h),2);
specs(end+1) = cell2struct({'Apass2','udouble',1,[],'magspec'},specfields(h),2);

% [EOF]
