%% Global estimation for data anaysis in thermometry
% 
% This code simulates the outcomes of mu energy measurements on a gas of n
% non-interacting spin-1/2 particles in thermal equilibrium. Then it processes 
% such outcomes to render a temperature estimate by using the theory of global 
% temperature estimation as put forward in J. Rubio et al (2020), arXiv:2011.13018.
%
% Running the code generates a plot with the result. 
% 
% Notes:
%   - dataOpt selects either (1) a fresh simulation, or (2) the simulated
%   dataset employed in Fig. 1b of arXiv:2011.13018. 
%   - In order to compare with the local methodology, this code also calculates
%   the temperature one would get from an optimal local estimator with
%   initial 'hint' at the temperature Ttrue-1. 
%   - The units have been chosen such that T is dimensionless.
%
% Jesús Rubio, PhD
% University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk
% Created: May 2021
% Last modified: June 2021

%% Initialisation
clear all

% Dataset option
dataOpt=1;

if dataOpt==1
    % New data set
    n=250; % number of particles
    mu=5000; % number of trials
    
elseif dataOpt==2
    % Data set in manuscript
    n=150;
    mu=500;
end

% Outcomes space
r=0:n;

%% Temperature estimation from simulated measurements: Bayesian approach

% Prior
Tmin=0.1; % lower limit
Tmax=10; % upper limit
dT=10^(-3);
dimT=round((Tmax-Tmin)/dT); % differential
T=linspace(Tmin,Tmax,dimT); % dimensionless temperature space
measure=sparse(1./T); % measure representing complete ignorance for scale parameters
prior=measure/trapz(T,measure); 

% Likelihood function
likelihood=zeros(n+1,length(T));
for xAux=1:n+1
    likelihood(xAux,:)=sparse(exp(-r(xAux)./T-n*log(1+exp(-1./T))+nCkLogJesus(n,r(xAux))));
end

%% Simulation

% 'True' temperature 
if dataOpt==1
    index_real=3900;
elseif dataOpt==2
    index_real=3900;
end

% Energy measurements (total energy of all n particles in each trial)
if dataOpt==1
    prob_sim=likelihood(:,index_real);
    rng('shuffle') % seed for the random generator
    for runs=1:mu
        
        auxiliar=cumsum(prob_sim)-rand;
        
        for x=1:n+1
            if auxiliar(x)>0
                outcome_index=x;
                break
            end
        end
        
        outcomes(runs)=r(outcome_index); %#ok<SAGROW>
        outcomesIndex(runs)=outcome_index; %#ok<SAGROW>
        
    end
    
elseif dataOpt==2
    outcomes=load('data_sample.txt');
    outcomesIndex=zeros(1,mu);
    for runs=1:mu 
        indexAux=find(outcomes(runs)==r);
       outcomesIndex(runs)=indexAux;
    end
end
        
% Inference
prob_temp=prior;
optEst=zeros(1,mu);optErr=zeros(1,mu);
locEst=zeros(1,mu);locErr=zeros(1,mu);
for runs=1:mu
        
    % Likelihood, joint, evidence and posterior functions
    joint=prob_temp.*likelihood(outcomesIndex(runs),:); % joint probability
    evidence=trapz(T,joint); % normalisation of Bayes theorem
    
    if evidence>1e-16
        posterior=joint/evidence; % posterior probability
    else
        posterior=0;
    end
    
    prob_temp=posterior; % this updates the posterior with the info of each new trial (Bayes theorem)
            
    % Optimal estimator
    aux=sparse(posterior.*log(T));
    optLogEst=trapz(T,aux);
    optEst(runs)=exp(optLogEst); % estimator in Eq.(6)
    
    % Optimal uncertainty
    optErr(runs)=trapz(T,aux.*log(T))-optLogEst^2;
    
    % Local theory
    localPoint=T(index_real)-1;
    locEst(runs)=localPoint+4*localPoint^2*mean(outcomes(1:runs))*cosh(1/(2*localPoint))^2/n-localPoint^2*(1+exp(-1/localPoint));   
    locErr(runs)=2*localPoint^2*cosh(1/(2*localPoint))/sqrt(n*runs);
    
end
optErrBar=optEst.*sqrt(optErr);

% Plots
shadedErrorBar(1:mu,optEst,optErrBar,'lineProps','b');
hold on
shadedErrorBar(1:mu,locEst,locErr,'lineProps','k');
plot(1:mu,T(index_real)*ones(1,mu),'r-','LineWidth',1.5)
hold off
fontsize=25;
xlabel('$\mu$','Interpreter','latex','FontSize',fontsize);
ylabel('$k_B\tilde{\theta}(\textbf{\emph{r}})/(\hbar \omega)$','Interpreter','latex','FontSize',fontsize);
legend('Opt. estimator, $\tilde{\vartheta}(\textbf{\emph{r}})$','Loc. estimator, $\tilde{\theta}_L(\textbf{\emph{r}})$','True temperature, $T$','Interpreter','latex','Location','southwest')
xlim([1 mu])
ylim([T(index_real)-1 T(index_real)+1])
set(gca,'FontSize',fontsize,'FontName','Times')
box on
grid
