function hblk = repeatingsequencestair(hTar,name)
%REPEATINGSEQUENCESTAIR


%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/11 15:49:57 $

error(nargchk(2,2,nargin,'struct'));

% check if simulink/Sources lib is available, if not load it.
issrclibloaded = 0;
srclibname = 'simulink';
srcblks_avail = issimulinkinstalled;
if srcblks_avail,
    wdsrcblk = warning;
    warning('off');
    if isempty(find_system(0,'Name',srclibname))
        issrclibloaded = 1;
        load_system(srclibname);
    end
    warning(wdsrcblk);
end

bname = 'simulink/Sources/Repeating Sequence Stair';
hblk = add_block(bname, [hTar.system '/' name]);

if issrclibloaded
    close_system(srclibname);
end



% [EOF]
