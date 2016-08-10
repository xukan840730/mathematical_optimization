function plugins = vhdlplugin
%VHDLPLUGIN   Plug-in file for the VHDL Filter product.

%   Author(s): J. Schickler
%   Copyright 2003-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.10 $  $Date: 2009/07/03 14:35:43 $

plugins.fdatool = @vhdldlg_plugin;

% -------------------------------------------------------------------------
function vhdldlg_plugin(hFDA)

% Add the VHDL menu item.
addtargetmenu(hFDA,'Generate HDL ...',{@vhdldlg_cb,hFDA},'generatehdl');

% -------------------------------------------------------------------------
function vhdldlg_cb(hcbo, eventStruct, hFDA) %#ok
Hd = getfilter(hFDA);
h = getcomponent(hFDA, '-class', 'hdlgui.fdhdldlg');

if isempty(h)
    [cando, msg] = ishdlable(Hd);
    if cando
        h = hdlgui.fdhdldlg(Hd);
        addcomponent(hFDA, h);
    else
        senderror(hFDA, msg);
        return;
    end
elseif ishandle(h.hHdlDlg)
    h.hHdlDlg.show; % bring the dialog to front
else
    h.hHdlDlg =  DAStudio.Dialog(h.hHdl);
end
% Create a listener on FDATool's Filter.
l = [...
    handle.listener(hFDA, 'FilterUpdated', {@lclfilter_listener, h.hHdl, h.hHdlDlg}); ...
    handle.listener(hFDA, 'CloseDialog', {@close_listener, h.hHdlDlg});...
    ];
setappdata(h, 'fdatool_listener', l);

% -------------------------------------------------------------------------
function lclfilter_listener(hFDA, eventData, hHdl, hDlg) %#ok

Hd = getfilter(hFDA);
[cando, msg] = hHdl.setfilter(Hd);

if ~isempty(hDlg) && ishandle(hDlg)
    hDlg.refresh;
end

if ~cando
    senderror(hFDA, msg);
end

% -------------------------------------------------------------------------
function close_listener(hFDA, eventData, hDlg) %#ok

if ~isempty(hDlg) && ishandle(hDlg)
    delete(hDlg);
end
% [EOF]
