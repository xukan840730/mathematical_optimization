function test = test(this, test)
%TEST Tests the M-code by running it.

%   Author(s): J. Schickler
%   Copyright 1988-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2008/04/21 16:31:33 $

if nargin < 2, test = testinit(class(this)); end

try
    eval(this.string);
    test = qeverify(test, {true, true});
catch ME
    disp(sprintf('M-code errored out with : %s', ME.message));
    test = qeverify(test, {true, false});
end

% [EOF]
