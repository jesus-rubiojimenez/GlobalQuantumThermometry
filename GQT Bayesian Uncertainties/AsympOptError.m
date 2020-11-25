function [eCramerRao]=AsympOptError(n,Tmin,Tmax)
% Given
%
%   - n free fermions in thermal equilibrium, and
%   - a prior temperature range [Tmin, Tmax],
% 
%   AsympOptError(n,Tmin,Tmax)
%
% returns the Bayesian mean logarithmic error in the large-n limit, which
% tunrs out to be Cramér-Rao-like quantity. 
%
% Notes:
%
%   - The units have been chosen such that T is dimensionless.
%   - The discretisation of the continuous variable T has been chosen such 
%     that the calculation is accurate for the interval [0.1, 10].
%   - While the code calculates the asymptotic error for any value of n
%     for comparison purposes, it is only meaningful when n>>1. 
%
% Jesús Rubio, PhD
% University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk
% Created: Sep 2020
% Last modified: Nov 2020

%% Parameter space
dT=10^(-2)/2;
dimT=round((Tmax-Tmin)/dT);
T=linspace(Tmin,Tmax,dimT);

%% Asymptotic Bayesian error (see paper for justification of this formula)
measure=sparse(1./T); % Complete ignorance for scale parameters
eCramerRao=4*trapz(T,T.*cosh(1./(2*T)).^2)/(trapz(T,measure)*n);

end