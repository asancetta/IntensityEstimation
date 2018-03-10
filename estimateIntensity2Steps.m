function [beta, hParms, gT] = ...
             estimateIntensity2Steps(XX,XatJump,jumpTime, constraintA, gStart)
% XatJump is matrix of covariates (number of jumps by number of covariates)
% evaluated at jump times only; XX is XDur'*X where X is the matrix of
% covariates observed over the whole sample (sample size by number of
% covariates) and XDur is X.*repmat(durations,1,K), where K is the number of columns in X;
% jumpTime is the time of a jump; 
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
ub            = 100;

%% --------- first step: baseline intensity estimation

durations      = [jumpTime(1);diff(jumpTime)];
% par=wblfit(durations);
% betaW=par(2);


[a,betaW]  = get_weibullParmsHt(durations,gStart);
hParms = betaW; 
hT                    = weibullMeasure(jumpTime,...
                                       ones(size(jumpTime))>0,...
                                       betaW,1);
  %  get_weibullLogLig(durations,1,1,gStart)
%% ---second step: etimation of g
A             = constraintA;
c             = zeros(size(A,1),1);

XY                    = XatJump'*(1./hT);

beta                 = quadprog(XX,-XY,A,c,[],[],...
                                0*ones(size(XY)),...
                                ub*ones(size(XY)));
beta(abs(beta)<smallNumber)  = 0;
gT                    = XatJump*beta;


end


function [a, b] = get_weibullParmsHt(x,g)

func       = @(parms) -1*get_weibullLogLig(x,parms(1),parms(2),g);
parmStart  = [1/mean(x),1];
A          = [];
b          = [];
Aeq        = [];
beq        = [];
lb         = [eps,eps];
ub         = [inf,inf];
parms      = fmincon(func,parmStart,A,b,Aeq,beq,lb,ub);

a          = parms(1);
b          = parms(2);


end

function out = get_weibullLogLig(x,a,b,g)

% g is a possibly time varing scale parameter in the Weibull with CDF(x) =
% 1-exp(-a*g.*x.^b) where g and x can have same dimnesion otherwise g is 1;

if nargin <4
   g  = 1;
end

out  = mean(log(a*b*g))+(b-1)*mean(log(x))-mean(a*g.*x.^b);

end

