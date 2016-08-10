function hblk = delay(hTar, name, latency, render)
%DELAY Add a Delay block to the model.
%   HBLK = DELAY(HTAR, NAME, LATENCY) adds a sum block named NAME, sets its
%   latency to LATENCY and returns a handle HBLK to the block.

%   Copyright 1995-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/11 15:49:23 $

error(nargchk(3,4,nargin,'struct'));
sys = hTar.system;

libname = 'simulink';
blockname = '/Discrete/Integer Delay';
latprop = 'NumDelays';
    
isloaded = 0;
w=warning;
warning('off');

if isempty(find_system(0,'Name', libname)),
    isloaded = 1;
    load_system(libname);
end
fullname = [libname blockname];

if nargin<4
    render=true;
end

if render
    hblk = add_block(fullname , [hTar.system '/' name], latprop, latency);
else         % then find the block and just update the block's parameters
    hblk1=find_system(sys,'SearchDepth',1,'Name',name);
    hblk=hblk1{1};
end

if isloaded,
    close_system(libname);
end

warning(w);
