function Hd = thisdesign(d)
%THISDESIGN  Local design method.

%   Author(s): R. Losada
%   Copyright 1988-2009 The MathWorks, Inc.
%   $Revision: 1.2.4.3.54.1 $  $Date: 2009/12/17 13:59:22 $


% Frequencies Have been prenormalized (0 to 1)

Hd = design(get(d, 'ResponseTypeSpecs'), d);
    
% Frequencies will be reset to what they were
