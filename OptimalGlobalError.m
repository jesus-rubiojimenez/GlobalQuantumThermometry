function [eOpt]=OptimalGlobalError(n,Tmin,Tmax)
% Given
%
%   - n non-interacting spin-1/2 particles in thermal equilibrium, and
%   - a prior temperature range [Tmin, Tmax],
% 
% the function
%
%     OptimalGlobalError(n,Tmin,Tmax)
%
% returns the minimium mean logarithmic uncertainty, as indicated by Eq.(7) 
% of J. Rubio et al. (2020), arXiv:2011.13018. Note that the chosen prior 
% represents complete ignorance within the interval [Tmin, Tmax]. 
%
% Remarks:
%
%   - The units have been chosen such that T is dimensionless.
%   - The discretisation of the continuous variable T has been chosen such 
%     that the calculation is accurate for the interval [0.1, 10].
%   - This code relies on the auxiliary function: nCkLog(n,k).
%
% Jesús Rubio, PhD
% University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk
% Created: Sep 2020
% Last modified: June 2021

%% Initialisation

if n>10^5
    warning('This code may not be accurate for systems with more than 10^5 particles.')
end

% Number of spin-1/2 particles
r=0:n;

% Parameter space
if Tmax>10
    warning('This code may be very slow when Tmax>>10.')
end
dT=10^(-3)/2;
dimT=round((Tmax-Tmin)/dT);
T=linspace(Tmin,Tmax,dimT);

%% Inference

% Prior
measure=sparse(1./T); % complete ignorance for a scale parameter   
prior=measure/trapz(T,measure); 

logInfo=0;
for index=1:length(r)
       
    % Likelihood, joint, evidence and posterior functions
    fermionLikelihood=sparse(exp(-r(index)./T-n*log(1+exp(-1./T))+nCkLog(n,r(index)))); 
    joint=prior.*fermionLikelihood;    
    evidence=trapz(T,joint);

    if evidence>1e-16
        posterior=joint/evidence;
    else
        posterior=0;
    end
        
    % Optimal logarithmic estimator
    aux1=sparse(posterior.*log(T));
    optLogEst=trapz(T,aux1);
    
    % Auxiliary term to calculate the mean logarithmic error 
    logInfo=logInfo+evidence*optLogEst^2;
end

% Optimal global mean logarithmic error
aux2=prior.*log(T).*log(T);
eOpt=trapz(T,aux2)-logInfo;
end