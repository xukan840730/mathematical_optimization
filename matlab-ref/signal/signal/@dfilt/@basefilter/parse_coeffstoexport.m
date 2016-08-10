function [hTar,domapcoeffstoports] = parse_coeffstoexport(Hd,hTar)
%PARSE_COEFFSTOEXPORT Store coefficient names and values into hTar for
%export.

%   Copyright 2009 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2009/07/14 04:00:50 $

state = hTar.MapCoeffsToPorts;

if strcmpi(state,'on')
    [mapstate coeffnames var] = mapcoeffstoports(Hd,'MapCoeffsToPorts','on',...
                                        'CoeffNames',hTar.CoeffNames);
                                
     % Coefficient names and variables 
     hTar.CoeffNames = coeffnames;
     setprivcoefficients(hTar,var);
end

domapcoeffstoports = strcmpi(state,'on');


% [EOF]
