function [s,str] = createStruct(h)
%CREATESTRUCT Create a structure for each type.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 2002/04/15 00:10:02 $

% Return a string for the name of the enumerated data type that will be created
str = 'iirgrpdelayTypes';

s(1).construct = 'filtdes.arbgrpdelay';
s(1).tag = 'Arbitrary group delay';
