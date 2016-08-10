function [errstr,P,f] = svwelch(x,Fs,valueArray)
%SVWELCH spectview Wrapper for Welch's method.
%  [errstr,P,f] = svwelch(x,Fs,valueArray) computes the power spectrum P
%  at frequencies f using the parameters passed in via valueArray:
%
%   valueArray entry     Description
%    ------------         ----------
%          1                Nfft
%          2                Window length
%          3                Window - integer
%                           index into        
%                   {'bartlett' 'blackman' 'boxcar' 'chebwin' ...
%                       'hamming' 'hanning' 'kaiser' 'triang'}
%          4                Window parameter (for chebwin and kaiser)
%          5                overlap - integer
%

%   Copyright 1988-2008 The MathWorks, Inc.
% $Revision: 1.12.4.1 $

errstr = '';
P = [];
f = [];

windowList = {'bartlett' 'blackman' 'boxcar' 'chebwin' ...
              'hamming' 'hanning' 'kaiser' 'triang'};

switch valueArray{3}
case {4,7}
    windStr = ...
    'window = feval(windowList{valueArray{3}},valueArray{2},valueArray{4});';
otherwise
    windStr = 'window = feval(windowList{valueArray{3}},valueArray{2});';
end

try
    eval(windStr);
catch ME
   errstr = {'Sorry, couldn''t evaluate window function; error message:'
      ME.message };
   return
end
   
nfft = valueArray{1};
noverlap = valueArray{5};

evalStr = '[P,f] = pwelch(x,window,noverlap,nfft,Fs);';

try
    eval(evalStr);
catch ME
    errstr = {'Sorry, couldn''t evaluate pwelch; error message:'
               ME.message };
    return
end

[P,f] = svextrap(P,f,nfft);
