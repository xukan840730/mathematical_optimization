function updateylabel(this, eventData)
%UPDATEYLABEL   Update the ylabel based on choice of frequency axis units.

%   Author(s): P. Pacheco
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2004/04/13 00:29:57 $

hprm = getparameter(this, getmagdisplaytag(this));
if isempty(hprm), return; end

ylabels = getylabels(this, eventData);

setvalidvalues(hprm, ylabels);

% if isempty(eventData) | strcmpi(eventData.Type,'NewValue'),
% 
%     %  Update the Ylabel if choice is density (i.e., /freq) to make sure
%     %  the ylabel is consistent with the frequency axis units.
% 
%     % Get the handle to the parameter object and get the Ylabel valid values.
%     hprm_magdisp  = getparameter(this, getmagdisplaytag(this));
%     Ylabel        = this.MagnitudeDisplay;
%     ylabelChoices = get(hprm_magdisp,'ValidValues');
%     
%     % Update the ylabel based on the new value of the x-axis units.
%     if ~isempty(ylabelChoices),  % Avoid initialization emptys
%         
%         hprm_freqmode = getparameter(this, 'freqmode');
%         normalizedmode = getsettings(hprm_freqmode, eventData);
%         
%         if strcmpi(normalizedmode,'off') & strcmpi(ylabelChoices{3},Ylabel),  % psd/rad/sample
%             % Change Ylabel to psd/Hz.  Set parameter object directly.
%             setvalue(hprm_magdisp,ylabelChoices{2});
%             
%         elseif strcmpi(normalizedmode,'on') & any(strcmpi(ylabelChoices(1:2),Ylabel)), % psd/hz,
%             % Change Ylabel to psd/rad/sample.  Set parameter object directly.
%             setvalue(hprm_magdisp,ylabelChoices{3});  
%         end
%     end
% end

% [EOF]
