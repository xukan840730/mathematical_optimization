function existPlugin = registerplugins(hFDA,mfile,structname)
%REGISTERPLUGINS Register mechanism for plug-ins.
%   REGISTERFDAPLUGINS(HFDA,MFILE,STRUCTNAME) search for paths containing 
%   MFILE and creates a plug-in specific structure called STRUCTNAME, 
%   see example below.
%   
%   existPlugin = REGISTERFDAPLUGINS(HFDA,MFILE,STRUCTNAME) stores a cell
%   array of structures with information about the plug-ins in object HFDA,
%   and returns a boolean flag indicating whether or not plug-ins exist.
%
%   EXAMPLE:
%  % Search and save plug-ins, if they exist, in FDATool.
%  existPlugins = registerfdaplugins(hFDA,'fdaregister.m','fdapluginstruct');

%   Author(s): P. Pacheco, P. Costa
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.8.4.1 $  $Date: 2008/05/31 23:28:40 $

% Initialized output variables in case of early return.
pluginList={};

toolboxes = which(mfile, '-all');
existPlugin = 1;
if isempty(toolboxes), 
    existPlugin = 0;
    return, 
end

% Plug-in structure name
structname = [structname,'.'];

j = 1;

for i=1:size(toolboxes, 1),
    
    % Use fopen to avoid checking out the toolbox license keys
    % This code is not robust to changes in the variable name "fdapluginstruct"
    fcnStr = '';
    fid = fopen(char(toolboxes(i)),'r'); 
    saveFlag = 0;
        
    while 1,
        fcnLine = fgetl(fid);
        if ~isstr(fcnLine), break, end
        if ~isempty(fcnLine) & ~strcmp(fcnLine(1),'%') & findstr(fcnLine,structname), 
            saveFlag = 1; 
        end
        if saveFlag,
            dotsIndex = findstr(fcnLine,'...');
            if isempty(dotsIndex),
                fcnStr=[fcnStr fcnLine 13];
            else
                % We're removing the "..." and literally concatenating the 
                % two affected lines to make one line
                fcnLine = strrep(fcnLine,'...','');
                % Note: this will remove the "..."s inside any quoted strings
                fcnStr=[fcnStr fcnLine];
            end
        end
    end
    fclose(fid);
    
    if ~isempty(fcnStr),

        % Initialize the structure
        eval([structname(1:end-1) ' = '''';']); 
        
        % fcnStr contains a string, evaluate this string to
        % creat the *pluginstruct structure
        eval(fcnStr);
        pluginStruct = eval([structname(1:end-1)]);
        if isfield(pluginStruct,'supportedplatforms'),
            platform = computer;
            supported = strmatch(platform,pluginStruct.supportedplatforms);
        else
            supported = 1;
        end
        if supported,
            pluginList{j} = pluginStruct;
            j = j + 1;
        end
    else
        % If there's no "*pluginstruct", just ignore the entry and keep on going...
    end
    
end

if isempty(pluginList),
    existPlugin = 0;
else
    % Store the plug-ins in "ud.plugins"
    hFig = get(hFDA,'figureHandle');
    registerPlugins(hFig,pluginList);
end

%------------------------------------------------------------------- 
function registerPlugins(hFig,pluginList)

ud = get(hFig,'Userdata');
ud.plugins.list = pluginList;  % Cell array of structures containing info about the plug-in.
set(hFig,'Userdata',ud)

% [EOF]
