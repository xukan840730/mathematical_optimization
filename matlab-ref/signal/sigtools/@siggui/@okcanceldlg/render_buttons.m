function render_buttons(hDlg)
%RENDER_BUTTONS Render the Dialog buttons (OK/Cancel)

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/13 00:24:42 $

% This can be a private method

hFig = get(hDlg,'FigureHandle');
sz   = dialog_gui_sizes(hDlg);

fsz  = figuresize(hDlg);
bgc  = get(0,'defaultuicontrolbackgroundcolor');

enabState = get(hDlg, 'Enable');

ctrlStrs = {sprintf('Cancel'),'OK'};
numbtns = 2;
uiWidth = largestuiwidth(ctrlStrs,'Pushbutton')+10*sz.pixf;

spacing = sz.uuhs*2;

okXPos     = fsz(1)/2-uiWidth*numbtns/2 - spacing*(numbtns-1)/2;
cancelXPos = okXPos + uiWidth + spacing;
applyXPos  = cancelXPos + uiWidth + spacing;

buttonPos = sz.button;

buttonPos([1,3]) = [okXPos uiWidth];

% NOTE: The converttdlg_cbs function updates the figure's userdata
cbs = dialog_cbs(hDlg);

% Render the "OK" pushbutton 
h.ok = uicontrol(hFig,...
    'Style','Push',...
    'BackgroundColor',bgc,...
    'Position',buttonPos,...
    'Visible','On',...
    'Enable',enabState, ...
    'String',ctrlStrs{2}, ...
    'Tag','dialog_ok', ...
    'Callback',cbs.ok);

buttonPos(1) = cancelXPos;

% Render the "Cancel" pushbutton 
h.cancel = uicontrol(hFig,...
    'Style','Push',...
    'BackgroundColor',bgc,...
    'Position',buttonPos,...
    'Visible','On',...
    'Enable',enabState,...
    'String',ctrlStrs{1},...
    'Tag','dialog_cancel',...
    'Callback',cbs.cancel);

h.warn = [];

set(hDlg,'DialogHandles',h);

% [EOF]
