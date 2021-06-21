function [C] = nCkLog(n,k)
% Given the integers n and k, with k in [n, 0), 
%
%   nCkLog(n,k) returns the
%
% natural logarithm of the binomial coefficient C(n, k), i.e., log[C(n, k)].
%
% This avoids numerical problems arising from the very large quantities 
% involved in the calculations with exponential probabilities, as is the
% case in quantum thermometry. 
%
% Jesús Rubio, PhD
% University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk
% Created: June 2020
% Last modified: Nov 2020
aux=0;
for x=k+1:n
    aux=aux+log(x)-log(x-k);
end

C=aux;
end
