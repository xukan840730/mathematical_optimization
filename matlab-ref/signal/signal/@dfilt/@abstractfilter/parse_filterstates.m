function hTar = parse_filterstates(Hd,hTar)
%PARSE_FILTERSTATES Store filter states in hTar for realizemdl

%   Copyright 2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/07/14 04:00:28 $

% Extract current filter states
IC = getinitialconditions(Hd);

% If the MapStates stored in hTar is not 'on', set the initial condition to
% 0.
if ~strcmpi(hTar.MapStates,'on')
    IC = zeros(size(IC));
end

% Store the filter states
setprivstates(hTar,IC);


% [EOF]
