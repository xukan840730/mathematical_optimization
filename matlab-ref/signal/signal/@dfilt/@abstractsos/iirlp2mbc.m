function [Ht,anum,aden] = iirlp2mbc(Hd, varargin)
%IIRLP2MBC IIR Lowpass to complex multiband transformation

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.2.4.2 $  $Date: 2005/06/16 08:17:16 $

[Ht,anum,aden] = ciirxform(Hd, @zpklp2mbc, varargin{:});

