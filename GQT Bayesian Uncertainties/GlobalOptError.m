function [eOpt]=GlobalOptError(n,Tmin,Tmax)
% Given
%
%   - n free fermions in thermal equilibrium, and
%   - a prior temperature range [Tmin, Tmax],
% 
%   GlobalOptError(n,Tmin,Tmax)
%
% returns the minimium mean logarithmic uncertainty associated with the
% optimal estimate, for a single shot, within the framework of global 
% quantum thermometry. In addition, the prior probability reflects complete 
% ignorance within the interval [Tmin, Tmax]. 
%
% Notes:
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
% Last modified: Nov 2020

%% Initialisation

if n>10^5
    warning('This code may not be accurate for systems with more than 10^5 free fermions.')
end

% Number of free fermions
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
measure=sparse(1./T); % Complete ignorance for a scale parameter   
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
    
    % Bayesian logarithmic information term
    logInfo=logInfo+evidence*optLogEst^2;
end

% Optimal mean logarithmic error
aux2=prior.*log(T).*log(T);
eOpt=trapz(T,aux2)-logInfo;
end