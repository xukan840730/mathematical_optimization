function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:10:35 $

if isdb(d),
    ab = get(d, 'Apass');
else
    ab = get(d, 'Epass');
end

if isprop(d, 'Order')
    nbands = get(d, 'Order');
else
    nbands = 2;
end

fn = getnyquist(d);
mu = get(d, 'MagUnits');

if nbands == 1,
    cmd{1}.amplitude  = ab;
    cmd{1}.filtertype = 'highpass';
    cmd{1}.magfcn     = 'cpass';
    cmd{1}.magunits   = mu;
else
    cmd.amplitude  = ab;
    cmd.filtertype = 'bandpass';
    cmd.magfcn     = 'cpass';
    cmd.magunits   = mu;
    
    cmd = repmat({cmd}, 1, floor(nbands/2));
    
    if rem(nbands, 2),
        cmd{end+1} = cmd{end};
        cmd{end}.filtertype = 'highpass';
    end
end

if ~isprop(d, 'Order')
    cmd{1}.filtertype = 'bandpass';
end

% [EOF]