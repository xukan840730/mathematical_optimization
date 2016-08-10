function thisstaticresponse(this, hax)
%THISSTATICRESPONSE   

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:13:56 $

if this.NormalizedFrequency, str = '.5';
else,                        str = 'Fs/4'; end

staticrespengine('drawpassband',   hax, [0   .45], [.9 1.1]);
staticrespengine('drawtransition', hax, [.45 .55]);
staticrespengine('drawstopband',   hax, [.55 1]);
staticrespengine('drawfreqlabels', hax, .5, str);

% [EOF]
