function varargout = render(this, varargin)
%RENDER Render the object

%   Author(s): J. Schickler
%   Copyright 1988-2008 The MathWorks, Inc.
%   $Revision: 1.12.4.9 $  $Date: 2009/01/05 18:01:14 $

if ~isrendered(this),
    
    sz = gui_sizes(this);
    origFontSize = get(0, 'defaultuicontrolfontsize');
    set(0, 'defaultuicontrolfontsize', sz.fontsize');
    
    if ~ispc
        origFontName = get(0, 'defaultuicontrolfontname');
        set(0, 'defaultuicontrolfontname', 'arial');
    end

    % Add the rendered properties
    add_properties(this);
    
    % Install listeners to these properties
    install_listeners(this);
    
    try
        % Allow subclasses to do whatever rendering they need.
        if nargout,
            [varargout{1:nargout}] = thisrender(this, varargin{:});
        else
            thisrender(this, varargin{:});
        end
    catch ME
        unrender(this);
        throwAsCaller(ME);
    end
    
    % Send the sigguiRendering event
    send(this, 'sigguiRendering', handle.EventData(this, 'sigguiRendering'));
    
    install_figurelistener(this);
    if isempty(this.Container)
        setunits(this, 'normalized');
    end
    set(0, 'defaultuicontrolfontsize', origFontSize);
    if ~ispc
        set(0, 'defaultuicontrolfontname', origFontName);
    end
end

% -------------------------------------------------------------------------
function add_properties(this)

% g108250 causes us to make all the properties public.

% Add a visible property that is public
p(1) = schema.prop(this, 'Visible', 'on/off');

% Add an enable property that is public
p(2) = schema.prop(this, 'Enable', 'on/off');

% Add a FigureHandle property that is publicget only
p(3) = schema.prop(this, 'FigureHandle', 'mxArray');
set(p(3), 'Visible', 'Off', 'SetFunction', @set_figurehandle, ...
    'GetFunction', @get_figurehandle);
% p(3).AccessFlags.PublicSet = 'off';

p(4) = schema.prop(this, 'Parent', 'mxArray');

p(5) = schema.prop(this, 'Container', 'mxArray');

% Add a Handles property that is private
p(6) = schema.prop(this, 'Handles', 'MATLAB array');
% p(4).AccessFlags.PublicSet = 'off';
% p(4).AccessFlags.PublicGet = 'off';

% Add a WhenRenderedListeners property that is private
p(7) = schema.prop(this, 'BaseListeners', 'handle.listener vector');
% p(5).AccessFlags.PublicSet = 'off';
% p(5).AccessFlags.PublicGet = 'off';

p(8) = schema.prop(this, 'WhenRenderedListeners', 'MATLAB array');
% p(6).AccessFlags.PublicSet = 'off';
% p(6).AccessFlags.PublicGet = 'off';

p(9) = schema.prop(this, 'RenderedPropHandles', 'schema.prop vector');

p(10) = schema.prop(this, 'Layout', 'mxArray');

set(this, 'RenderedPropHandles', p);
set(this, 'Enable', 'On');
set(this, 'Parent', -1);

set(p(6:10), 'Visible', 'Off');

% -------------------------------------------------------------------------
function hf = set_figurehandle(this, hf)

set(this, 'Parent', hf);

hf = [];


% -------------------------------------------------------------------------
function hf = get_figurehandle(this, hf) %#ok

hf = ancestor(this.Parent, 'figure');
if isempty(hf)
    hf = -1;
end

% -------------------------------------------------------------------------
function install_listeners(this)

% Create listeners on the visible and enable properties
listener(1) = handle.listener(this, this.findprop('Visible'), ...
    'PropertyPostSet', @lclvisible_listener);
listener(2) = handle.listener(this, this.findprop('Enable'), ...
    'PropertyPostSet', @lclenable_listener);

% Set the callback target to the object
set(listener, 'CallbackTarget', this);

% Save the listeners in WhenRenderedListeners
set(this, 'BaseListeners', listener);


% -------------------------------------------------------------------------
function install_figurelistener(this)

hFig = get(this, 'Parent');

if ishghandle(hFig),

    listener = get(this, 'BaseListeners');
    
    % Create the listener
    listener(end+1) = handle.listener(hFig, 'ObjectBeingDestroyed', @lclfigure_listener);
    listener(end+1) = handle.listener(this, 'ObjectBeingDestroyed', @objectbeingdestroyed_listener);
    set(listener,'CallbackTarget',this);
    
    set(this, 'BaseListeners', listener);
end

% -------------------------------------------------------------------------
%       Local Listeners (local function handles are faster)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function lclfigure_listener(this, eventData)

% Make sure the object is still valid.  The event is being cached and the
% object is being deleted before coming into this listener
if isa(this, 'siggui.siggui'),
    figure_listener(this, eventData)
end

% -------------------------------------------------------------------------
function lclenable_listener(this, eventData)

enable_listener(this, eventData)

% -------------------------------------------------------------------------
function lclvisible_listener(this, eventData)

visible_listener(this, eventData)

% -------------------------------------------------------------------------
function objectbeingdestroyed_listener(this, eventData) %#ok

if isrendered(this) && ~strcmpi(this.Tag, 'siggui.cfi')
    set(this, 'Visible', 'Off');
    unrender(this);
end

% [EOF]
