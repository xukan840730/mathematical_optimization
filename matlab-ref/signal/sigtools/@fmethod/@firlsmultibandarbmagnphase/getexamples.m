function examples = getexamples(this)
%GETEXAMPLES   Get the examples.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2005/08/20 13:28:48 $

examples = {{...
    'Bandpass Filter with a Low Group Delay',...
    'F1 = linspace(0,.25,30);',...
    'F2 = linspace(.3,.56,40);',...
    'F3 = linspace(.62,1,30);',...
    'gd = 12;',...
    'H1 = zeros(size(F1));',...
    'H2 = exp(-j*pi*gd*F2);',...
    'H3 = zeros(size(F3));',...
    'f=fdesign.arbmagnphase(''N,B,F,H'',40,3,F1,H1,F2,H2,F3,H3);',...
    'Hd = design(f,''firls'');',...
    'fvtool(Hd,''DesignMask'', ''on'')',...
    },...
    {'Use Weights to Improve Stopband Characteristics',...
    'W1 = 10*ones(size(F1));',...
    'W2 = ones(size(F2));',...
    'W3 = 10*ones(size(F3));',...
    'Hd(2) = design(f,''firls'',''B1Weights'',W1,''B2Weights'',W2,''B3Weights'',W3);',...
    'fvtool(Hd,''DesignMask'', ''on'')',...
    }};


% [EOF]
