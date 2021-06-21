function [eCramerRao]=LocalLimitError(n,Tmin,Tmax)
% Given
%
%   - n non-interacting spin-1/2 particles in thermal equilibrium, and
%   - a prior temperature range [Tmin, Tmax],
% 
% the function
%
%   LocalLimitError(n,Tmin,Tmax)
%
% returns the mean logarithmic error in the limit of local prior information, 
% as indicated by Eq.(12) of J. Rubio et al. (2020), arXiv:2011.13018. This 
% turns out to be a Cramér-Rao-like quantity. 
%
% Notes:
%
%   - The units have been chosen such that T is dimensionless.
%   - The discretisation of the continuous variable T has been chosen such 
%     that the calculation is accurate for the interval [0.1, 10].
%   - The code calculates this error for any value of n in order to compare
%     with the true optimum. However, note that the specific value of this
%     local error is meaningful only when n >> 1. 
%
% Jesús Rubio, PhD
% University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk
% Created: Sep 2020
% Last modified: June 2021

%% Parameter space
dT=10^(-2)/2;
dimT=round((Tmax-Tmin)/dT);
T=linspace(Tmin,Tmax,dimT);

%% Local Bayesian error (see Eq.(12) of manuscript)
measure=sparse(1./T); % complete ignorance for scale parameters
eCramerRao=4*trapz(T,T.*cosh(1./(2*T)).^2)/(trapz(T,measure)*n);
end
