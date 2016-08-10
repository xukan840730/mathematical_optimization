function sethdl_abstractfilter(this, hhdl)
%SETHDLPROPSBASEFILTER Set the common props for HDLFILTER  from filter
%object

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2007/09/18 02:30:16 $

hhdl.FilterStructure = this.FilterStructure;

hhdl.InputSLType = conv2sltype(q, 'InputWordLength', 'InputFracLength', true);
hhdl.OutputSLType = conv2sltype(this.Filterquantizer, 'OutputWordLength', 'OutputFracLength', true);

% copy castbeforesum from filterobj to hdlfiltercomp if it is available in
% filterobj otherwise set to false (default)
% this is mostly used in case these objects are used in cascade stages
if any(strcmpi(fieldnames(get(this)),'Arithmetic')) && ...
        strcmpi(this.arithmetic,'fixed')
    if any(strcmpi(fieldnames(this),'CastBeforeSum'))
        if this.CastBeforeSum,
            hhdl.Castbeforesum = this.castbeforesum;
        else
            hhdl.Castbeforesum = false;
        end
    end
end

% [EOF]
