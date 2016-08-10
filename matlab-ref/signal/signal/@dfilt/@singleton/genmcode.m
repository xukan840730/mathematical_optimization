function str = genmcode(Hd, objname, place) %#ok
%GENMCODE Generate the M-code to reproduce the filter.
%   GENMCODE(Hd, NAME) Generate M-code to reproduce the filter which will
%   be stored in the variable specified by the string NAME.  If NAME is not
%   specified 'Hd' is used.

%   Author(s): J. Schickler
%   Copyright 1988-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2006/06/27 23:34:57 $

if nargin < 2, objname = 'Hd'; end

object = class(Hd);

labels = propnames(Hd);
coeffs = propvalues(reffilter(Hd));

for indx = 1:length(coeffs)
    coeffs{indx} = genmcodeutils('array2str', coeffs{indx}, 16);
    
    % Remove all the extra spaces.
    [s, f] = regexp(coeffs{indx}, ' + ');
    idx = [];
    for jndx = 1:length(s)
        idx = [idx s(jndx)+1:f(jndx)];
    end
    
    coeffs{indx}(idx) = [];
    descs{indx} = sprintf('%s coefficient vector', labels{indx});
    labels{indx}(strfind(labels{indx}, ' ')) = '_';
end

% Format the labels and coeffs.
inputs = sprintf('%s, ', labels{:});
inputs(end-1:end) = [];

h = sigcodegen.mcodebuffer;

h.add(h.formatparams(labels, coeffs, descs));
h.cr;
h.cradd('%s = %s(%s);', objname, object, inputs);

if isfdtbxinstalled && isprop(Hd, 'Arithmetic')
    
    h.cr;
    h.cradd('set(%s, ''Arithmetic'', ''%s''', objname, Hd.Arithmetic);

    hq = genmcode(Hd.filterquantizer);
    if ~isempty(hq)
        h.add(hq);
    end
    h.add(');');
end

str = h.string;

% [EOF]
