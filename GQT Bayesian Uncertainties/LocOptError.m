function [eLocOpt]=LocOptError(n,Tmin,Tmax)
% Given
%
%   - n free fermions in thermal equilibrium, and
%   - a prior temperature range [Tmin, Tmax],
% 
%   LocOptError(n,Tmin,Tmax)
%
% returns the Bayesian mean logarithmic uncertainty associated with the
% locally-optimal estimate proposed in Phys Rev E, 83 011109 (2011) for this
% fermionic system (see the documentation of LocOptEstimator(r,n,Tmin,Tmax)). 
% In addition, the prior probability reflects complete ignorance within the 
% interval [Tmin, Tmax]. 
%
% Notes:
%
%   - The units have been chosen such that T is dimensionless.
%   - The discretisation of the continuous variable T has been chosen such 
%     that the calculation is accurate for the interval [0.1, 10].
%   - This code relies on the auxiliary functions: nCkLog(n,k) and 
%     LocOptEstimator(r,n,Tmin,Tmax).
%
% Jesús Rubio, PhD
% University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk
% Created: Sep 2020
% Last modified: Nov 2020

%% Initialisation

if n>10^5
    warning('This code may not be accurate for systems with more than 10^5 free spins.')
end

% Number of free fermions
r=0:n;

% Parameter space
dT=10^(-3)/2;
dimT=round((Tmax-Tmin)/dT);
T=linspace(Tmin,Tmax,dimT);

%% Inference

eMLEenergyLoc=0;
for indexr=1:length(r)
    
    % Likelihood
    fermionLikelihood=exp(-r(indexr)./T-n*log(1+exp(-1./T))+nCkLog(n,r(indexr)));
    
    % Deviation function for given r and T
    errFun=(log(LocOptEstimator(r(indexr),n,Tmin,Tmax))-log(T)).^2;
    
    % Mean logarithmic error (temperature-dependent; local version)
    eMLEenergyLoc=eMLEenergyLoc+fermionLikelihood.*errFun;
end

% Prior
measure=sparse(1./T); % Complete ignorance for a scale parameter 
prior=measure/trapz(T,measure);

% Mean logarithmic error (averaged over T; global version)
eLocOpt=trapz(T,prior.*eMLEenergyLoc);
end