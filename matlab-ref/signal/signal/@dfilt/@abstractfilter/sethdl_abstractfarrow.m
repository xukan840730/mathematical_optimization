function sethdl_abstractfarrow(this, hhdl)
%SETHDL_ABSTRACTFARROW Set the properties of hdlfiltercomp (hhdl) from the
%filter object.

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2007/09/18 02:30:07 $
this.sethdl_abstractfilter(hhdl);
[hhdl.RoundMode, hhdl.OverflowMode] = conv2hdlroundoverflow(this);

hhdl.ProductSLType = conv2sltype(this.filterquantizer, 'ProductWordLength', 'ProductFracLength', true);
hhdl.AccumSLType = conv2sltype(this.filterquantizer, 'AccumWordLength', 'AccumFracLength', true);

hhdl.FDSLType = conv2sltype(this.filterquantizer, 'FDWordlength', 'FDFraclength', false);
% [EOF]
