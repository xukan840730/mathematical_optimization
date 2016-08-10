function pos = parserenderinputs(this, varargin)
%PARSERENDERINPUTS Parse for the inputs to render
%   PARSERENDERINPUTS Parse for the inputs to render (hFig, position)

%   Author(s): J. Schickler
%   Copyright 1988-2008 The MathWorks, Inc.
%   $Revision: 1.3.4.5 $  $Date: 2009/01/05 18:01:11 $

hFig = -1;
pos  = [];

for i = 1:length(varargin)
    if isnumeric(varargin{i}), 
        if ishghandle(varargin{i})
            
            if ishghandle(varargin{i}, 'figure') || ...
                    ishghandle(varargin{i}, 'uipanel') || ...
                    ishghandle(varargin{i}, 'uicontainer')
                hFig = varargin{i};
            end
        elseif length(varargin{i}) == 4,
            pos = varargin{i};
        end
    end
end

if ~ishghandle(hFig),
    if ishghandle(this.Parent, 'figure') || ...
            ishghandle(this.Parent, 'uipanel') || ...
            ishghandle(this.Parent, 'uicontainer')
        hFig = this.Parent;
    else
        hFig = gcf;
    end
end

set(this, 'Parent', hFig);

% [EOF]
