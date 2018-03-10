function [beta, hParms] = ...
             estimateIntensity2Steps(XX,XatJump,time, isJump, constraintA, gStart)
% XatJump is matrix of covariates (number of jumps by number of covariates)
% evaluated at jump times only; XX is XDur'*X where X is the matrix of
% covariates observed over the whole sample (sample size by number of
% covariates) and XDur is X.*repmat(durations,1,K), where K is the number of columns in X;
% time is the time of all updates including jump times;
% is jump is True if the corresponding row in time is a jump time, i.e. time(isJump) are the jump times 
% constraintA: constraint to use in the estimation, it is the output of the function 
         
         
         
if nargin < 4
    constraintA  = [];
end

if nargin < 5
    gStart  = 1;
end


%% hard coded parameters: change them is necessary

% set beta coefficients in g estimation equal to zero if smaller than this
smallNumber   =  1e-6;
% maximum value for beta coefficients:
ubBeta        = 100;

%% --------- first step: baseline intensity estimation
  
jumpTimes          = makeItLeftContinuous(time,time,isJump);
jumpTimes(1)       = 0; %assume there is a jump at time 0
timeLag            = [nan(1,1); time(1:end-1)];
R                  = (time-jumpTimes);  
Rlag               = (timeLag-jumpTimes);
Rlag(1)            = 0; 
RJump              = R(isJump);      
func  = @(x) -1*get_weibullIntensityLikelihood(x(1), x(2),...
                                                 R,Rlag,RJump, gStart);
parmStart = [1,1/sqrt(mean(gStart.^2))];
lb        = [.0001,.0001];
ub        = [10,100];
parm      = fmincon(func,parmStart,[],[],[],[],lb,ub);
                                             
betaW=parm(1);
hParms = betaW; 
flagDensity           = 1;% compute the hazrd and not the integral                                           
hT                    = weibullMeasure(time(isJump),ones(sum(isJump),1)>0,betaW,flagDensity);
  
%% ---second step: etimation of g
A             = constraintA;
c             = zeros(size(A,1),1);

XY                    = XatJump'*(1./hT);

beta                 = quadprog(XX,-XY,A,c,[],[],...
                                0*ones(size(XY)),...
                                ubBeta*ones(size(XY)));
beta(abs(beta)<smallNumber)  = 0;

end


function out = get_weibullIntensityLikelihood(shape, scale, R, Rlag,RJump, g)
% lambda = scale*shape*R^(shape-1)*g
out     = size(RJump,1)*log(scale*shape)+(shape-1)*sum(log(RJump))...
                          -scale*sum(((R.^shape)-(Rlag.^shape)).*g);

                      
end
