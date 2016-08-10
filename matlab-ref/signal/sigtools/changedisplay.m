function str = changedisplay(s, varargin)
%CHANGEDISPLAY   Change the display.
%   CHANGEDISPLAY(S, PROP, VAL) Change the display of structure S so that
%   field PROP shows the string VAL without the extra ''.
%
%   CHANGEDISPLAY(S, PROP1, VAL1, PROP2, VAL2, etc) Change the display for
%   multiple prop/val pairs.
%
%   s.field1 = true;
%   s.field2 = 10;
%
%   disp(s);
%
%   disp(changedisplay(s, 'field1', 'true'));

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2004/12/26 22:23:15 $

for indx = 1:2:length(varargin)
    s.(varargin{indx}) = varargin{indx+1};
end

% Convert the structure to a string so that we can manipulate it.
str = evalc('s');

for indx = 2:2:length(varargin)
    
    % Remove the extra '' around the string.
    str = strrep(str, sprintf('''%s''', varargin{indx}), varargin{indx});

    % Remove the first line since it reads 's = \n'
    sndx = strfind(str, 's =');
    str(1:sndx+2) = [];
end

% Look for extra new line feeds, we want to get rid of them all.
sndx = min(regexp(str, '\w'));
sndx = find(str(1:sndx-1) == char(10), 1, 'last');
str(1:sndx) = [];

sndx = find(str~=char(10), 1, 'last');
str(sndx+2:end) = [];

% [EOF]
