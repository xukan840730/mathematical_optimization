function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the designobj.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2005/10/14 16:28:29 $

if isfdtbxinstalled
    
    %#function fdfmethod.kaiserhbastop
    %#function fdfmethod.eqriphbastop
    %#function fdfmethod.elliphalfbandastop
    %#function fdfmethod.iirhalfbandeqripastop
    designobj.kaiserwin  = 'fdfmethod.kaiserhbastop';
    designobj.equiripple = 'fdfmethod.eqriphbastop';
    designobj.ellip      = 'fdfmethod.elliphalfbandastop';
    designobj.iirlinphase = 'fdfmethod.iirhalfbandeqripastop';
else
    designobj = [];
end

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
