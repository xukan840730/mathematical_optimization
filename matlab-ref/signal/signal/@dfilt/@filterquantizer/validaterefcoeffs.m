function validaterefcoeffs(q, prop, val)
%VALIDATEREFCOEFFS   

%   Author(s): V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2004/04/12 23:58:01 $

if ~(strcmpi(class(val), 'double') || strcmpi(class(val), 'single') ...
        || strcmpi(class(val), 'embedded.fi') ...
        || strncmpi(class(val), 'int', 3) || strncmpi(class(val), 'uint', 4)),
    error(generatemsgid('invalidDataType'), [prop ' must be of class fi, double, int* or uint*.']);
end

if issparse(val),
    error(generatemsgid('invalidDataType'), [prop ' cannot be a sparse matrix.']);
end

% [EOF]
