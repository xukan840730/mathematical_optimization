function genmcode(this, file)
%GENMCODE Generate M-code

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.4.5 $  $Date: 2004/12/26 22:23:05 $

mcode = get(this, 'MCode');

% If there is no mcode, generate the default mcode.
if isempty(ishandle(mcode)),
    [filt, mcode] = defaultfilter(this);
end
if isempty(mcode),
    error(generatemsgid('noMcode'), 'There is no M code to generate a file.');
end

if nargin < 2,
    [file, path] = uiputfile('*.m', 'Generate M-file', 'untitled.m');
    if isequal(file, 0),
        return;
    end
    file = fullfile(path, file);
    if isempty(strfind(file, '.')),
        file = [file '.m'];
    end
end

% Set up the options for the public writer.
if isa(getfilter(this), 'mfilt.abstractmultirate')
    opts.H1         = 'Returns a multirate filter object.';
else
    opts.H1         = 'Returns a discrete-time filter object.';
end
opts.outputargs = 'Hd';

% Call the public writer.
genmcode(file, mcode, opts);

% Launch the editor with the file.
edit(file);

% [EOF]
