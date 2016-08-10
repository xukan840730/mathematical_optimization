function specs = fp1_whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 2002/04/15 00:38:15 $

% Prop name, data type, default value, listener callback
specs = cell2struct({'Fpass','udouble',9600,[],'freqspec'},specfields(h),2);



