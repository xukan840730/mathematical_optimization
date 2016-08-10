function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the designobj.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2007/10/23 18:49:22 $

if isfdtbxinstalled
    
    %#function fmethod.elliplpfstop
    %#function fdfmethod.eqriplpfstop
    designobj.ellip      = 'fmethod.elliplpfstop';
    designobj.equiripple = 'fdfmethod.eqriplpfstop';
else
    designobj = [];
end

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
