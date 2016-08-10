function [s,c] = bintheor(x,y,n)
%BINTHEOR Binomial theorem.
% This m-file gives the expansion of powers of sums of any real or complex
% numbers x and y, and any nonnegative integer n. It is also known as the
% Newton's binomial. From it emerges the discrete binomial (positive)
% distribution. The formulation is as,
%                    
%                   n_
%                   \       n!
%     (x +/- y)^n = /_   -------- x^(n-k) y^k
%                   k=0  k!(n-k)!
%
%
% The
%         n!                                        
%      -------- , is the called binomial coefficient. 
%      k!(n-k)!
%
% Syntax: function bintheor(x,y,n) 
%      
%     Input:
%          x,y - pair of interested terms to expand
%            n - coefficient/power to increase the binomial theorem
%     Output:
%          result of the Binomial theorem sum (default)
%          vector of the binomial theorem values (optional)
%
% Example: For x=3, y=2, n=7
%
% Calling on Matlab the function: 
%             [s,c]=bintheor(3,2,7)
%
%  Answer is:
%
%  s = 7.8125e+004
%
%  c =
%  1.0e+004 *
%
%  0.2187    1.0206    2.0412    2.2680    1.5120    0.6048    0.1344    0.0128 
%
% Created by A. Trujillo-Ortiz, R. Hernandez-Walls, K. Barba-Rojo
%            and N. Nunez-Valencia
%            Facultad de Ciencias Marinas
%            Universidad Autonoma de Baja California
%            Apdo. Postal 453
%            Ensenada, Baja California
%            Mexico.
%            atrujo@uabc.mx
% Copyright. November 12, 2008.
%
%  To cite this file, this would be an appropriate format:
%  Trujillo-Ortiz, A., R. Hernandez-Walls, K. Barba-Rojo and 
%    N. Nunez-Valencia. (2008). bintheor:Binomial theorem. A MATLAB file.
%    [WWW document]. URL http://www.mathworks.com/matlabcentral/fileexchange/
%    loadFile.do?objectId=22085
%
% Reference:
% M. Abramowitz and Stegun, I.A. (1972), Handbook of Mathematical
%      Functions. Government Printing Office, p.10.
%

if  nargin < 3
    error('TooFewInputs:BINTHEOR requires three input arguments.');
end

c = [];
for i =0:n,
    ni = gammaln(n + 1) - gammaln(i + 1) - gammaln(n - i + 1);
    lni = ni + (n - i).*log(x) + i.*log(y);
    yi = exp(lni);
    c = [c yi];
    c = real(c);
end

s = sum(c);

return