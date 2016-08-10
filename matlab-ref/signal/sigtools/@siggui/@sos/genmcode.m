function str = genmcode(hObj)
%GENMCODE Generate M-code

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2003/03/02 10:28:37 $

% Map proper SOS scaling strings 
filtobj = get(hObj, 'Filter');
values = {sprintf('''%s''', get(hObj, 'Direction'))};

if isa(filtobj, 'dfilt.df2') | isa(filtobj, 'dfilt.df2sos'),
    
    params = {'scale', 'dir'};
    descs  = {'Scaling', 'Direction Flag'};
    
    scale = get(hObj, 'Scale');
    inputs = ', scale';
    switch scale
        case 'L-2',
            values = {'2', values{:}};
        case 'L-infinity',
            values = {'inf', values{:}};
        otherwise
            values = {'0', values{:}};
    end
else
    params = {'dir'};
    descs  = {'Direction Flag'};
    inputs = '';
end

if isa(filtobj, 'dfilt.abstractsos'),
    comments = '% Scale the second-order sections filter.';
else
    comments = '% Convert the filter to second-order sections.';
end

str = { ...
        genmcodeutils('formatparams', params, values, descs), ...
                '', ...
                comments, ...
        sprintf('Hd = sos(Hd, dir%s);', inputs), ...
    };

str = sprintf('%s\n', str{:}); str(end) = [];

% [EOF]
