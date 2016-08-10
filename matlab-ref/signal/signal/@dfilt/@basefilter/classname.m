function fstruct = classname(Hb)
%CLASSNAME Returns the class name of the object
%   CLASSNAME(Hb) Returns the class name of the object (the abbreviated 
%   filter structure)

% Author(s): J. Schickler
% Copyright 1988-2002 The MathWorks, Inc.
% $Revision: 1.1 $ $Date: 2002/07/18 21:44:12 $

fstruct = get(Hb.classhandle, 'Name');

% [EOF]
