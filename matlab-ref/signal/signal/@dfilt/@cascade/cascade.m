function Hd = cascade(varargin)
%CASCADE Create a cascade of discrete-time filters.
%   Hd = CASCADE(Hd1, Hd2, etc) returns a discrete-time filter, Hd, of type
%   cascade, which is a serial interconnection of two or more dfilt filters,
%   Hd1, Hd2, and so on. For more information on filter objects, see the 
%   <a href="matlab:web([matlabroot,'\toolbox\signal\sigdemos\html\dfiltdemo.html'])">Getting Started with Discrete-Time Filters</a> demo. The block diagram of
%   this cascade looks like:
%
%      x ---> Hd1 ---> Hd2 ---> etc. ---> y
%
%   Note that with the Filter Design Toolbox installed, one usually does
%   not construct CASCADE filters explicitly. Instead, one obtains these
%   filters as a result from a design using <a href="matlab:help fdesign">FDESIGN</a>. 
%
%   % EXAMPLE #1: Direct instantiation
%   Hd = dfilt.dffir([0.05 0.9 0.05]);
%   Hgain = dfilt.scalar(2);
%   Hcas = dfilt.cascade(Hgain,Hd)
%   realizemdl(Hcas)    % Requires Simulink
%   
%   % EXAMPLE #2: Design an Interpolated FIR lowpass filter 
%   Hcas = design(fdesign.lowpass('Fp,Fst,Ap,Ast',.1, .12, 1, 60), 'ifir')
%   fvtool(Hcas)        % Analyze filter
%   x = randn(100,1);   % Input signal
%   y = filter(Hcas,x); % Apply filter to input signal
%
%   See also DFILT/STRUCTURES   
  
%   Copyright 1988-2009 The MathWorks, Inc.
%   $Revision: 1.9.4.9 $  $Date: 2009/11/13 05:03:21 $

if nargin == 0,
    varargin = {dfilt.dffir(1),dfilt.dffir(1)};
end

Hd = dfilt.cascade;

Hd.FilterStructure = 'Cascade';

% Check that all are dfilts before starting to set parameters.
for k=1:length(varargin)
  if isnumeric(varargin{k})
    g = squeeze(varargin{k});
    if isempty(g) || length(g)>1
      error(generatemsgid('Empty'),'Numeric section in a cascade must be a nonempty scalar.');
    end
    varargin{k} = dfilt.scalar(g);
  end
  if isfdtbxinstalled,
      if ~(isa(varargin{k}(end),'dfilt.singleton') || ...
              isa(varargin{k}(end),'dfilt.multistage') || ...
              isa(varargin{k}(end),'dfilt.abstractfarrowfd'))
          error(generatemsgid('DFILTErr'),'Cascades must be made of discrete-time (DFILT).');
      end
  else
      if ~(isa(varargin{k}(end),'dfilt.singleton') || isa(varargin{k}(end),'dfilt.multistage'))
          error(generatemsgid('DFILTErr'),'Cascades must be made of discrete-time (DFILT) objects.');
      end
  end
end

for k=1:length(varargin)
  Hd.Stage = [Hd.Stage; varargin{k}(:)];
end
