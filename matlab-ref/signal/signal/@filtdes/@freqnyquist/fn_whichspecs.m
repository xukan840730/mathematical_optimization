function specs = fn_whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.2.4.2 $  $Date: 2005/06/16 08:26:04 $



% Prop name, data type, default value, listener callback
specs(1) = cell2struct({'Band','spt_uint32',4,[],''},specfields(h),2);
