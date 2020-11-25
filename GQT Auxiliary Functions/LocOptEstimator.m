function [TLocOptEst]=LocOptEstimator(r,n,Tmin,Tmax)
% Given
%
%   - n free fermions in thermal equilibrium, and
%   - a prior temperature range [Tmin, Tmax],
%
% suppose we measure the number of fermions that are in their excited state, 
% finding the outcome r. Then,
% 
%   LocOptEstimator(r,n,Tmin,Tmax)
%
% maps the outcome r to and returns a temperature estimate that is optimal
% according to local thermometry (but not within the framework of global
% estimation theory). In particular, this estimator  arises 
% from inverting the functional relation between average energy and 
% temperature in statistical mechanics, as prescribed in Phys Rev E, 83 
% 011109 (2011).
%
% Notes:
%
%   - The units have been chosen such that T is dimensionless.
%   - The truncation for the estimation to lie in the interval [Tmin, Tmax] 
%     is needed when using a Bayesian approach. 
%
% Jesús Rubio, PhD
% University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk
% Created: Sep 2020
% Last modified: Nov 2020

%% Locally-optimal estimator for a given measurement outcome r
if n/(exp(1/Tmax)+1)>r && r>0
    TLocOptEst=1/log(n/r-1);
elseif r==0 || r>=n/2
    TLocOptEst=Tmin;
elseif (n/2>r && r>n/(exp(1/Tmax)+1)) || n==2*r
    TLocOptEst=Tmax;
end

end