%Example of work flow

%% Generate a sample from the intensity .1+x1.^a+x2.^b to show how to use the code 
K            = 2; % number of covariates
constPoly    = 0.01;
expPoly      = [1,3];
gFunc        = @(x) constPoly+x(:,1).^expPoly(1)+x(:,2).^expPoly(2); %the intensity function
%gFunc        = @(x) .1+sum(x.^2,2); %another the intensity function
n            = 1e6;% sample size
X            = rand(n,K);% the covariates
g            = gFunc(X); %intensity evaluated at the covariates
durations    = -log(rand(n,1))./g;
jumpTimes    = cumsum(durations);
%all times are jump times; generally, time includes time of updates in X's,
% not just jump times; 
time         = jumpTimes; %all times are jump times in this example
isJump       = ones(size(time))>0; %all times are jump times

%% NOTE: TIME MUST NOT START AT ZERO OR SOME ARBITRARY TIME
%% IT NEEDS TO BE THE CUMSUM OF DURATIONS AS DONE ABOVE

%% prepare the data for estimation

% set the parameters for estimation and define the Bernstein basis 
orderPoly  = 8;%order of Bernstein polynomail
basis      = get_bernsteinBasisMultivariate(orderPoly,X);   

%isJump: 1 if the row is a jump, zero otherwise; 
XX        = bsxfun(@times,basis,durations)'*basis;
XatJump   = basis(isJump,:);

%identify the covariate in basis: used to construct a constraint for a
%specific covariate
indVar     = kron((1:K),ones(1,orderPoly+1));

%% estimate the parameters

% define the constraint, e.g. convex monotone increasing for variable 2

k      = 2;%Variable 2 only
a      = [1,-1];% for monotone increasing;
%a      = -[1,-1];% for monotone decreasing;
%a      = [-1,2,-1];% for convex;
%a      = -[-1,2,-1];% for concave;
[A1, c1] = get_constr(indVar,k,a);% check function for details

a        = [-1,2,-1];% for convex;
[A2, c2] = get_constr(indVar,k,a);
A        = [A1;A2];%put constraints together: monotone increasing and convex 
c        = [c1;c2]; % this is always zero Ab<=c and c=0 (vector)

% carry out the estimation
gT  = 1;% start with guess of gT=1
s  =0;
maxIter  = 5;% iterate 5 times

while s<maxIter
    
%     [beta, hParm, gT] = ...
%              estimateIntensity2Steps_test(XX,XatJump,jumpTimes, A,gT/mean(gT));
% 
    [beta, hParm] = ...
             estimateIntensity2Steps_test(XX,XatJump,time, isJump, A,gT/mean(gT));
    gT                = basis*beta;         
    s =s+1;         
end
%% compute the integrated intensity 
[Lambda, gT ]  = get_integratedIntensity(basis, ...
                                        time,...
                                        isJump,...
                                        hParm,...
                                        beta);

%% plot impact function for each variable
close all;         
figure(1);
varId    = 1;%plot the intensity component for the first variable
bernsteinPolyPlot(beta,orderPoly,varId);%plot the estimated function

hold on; % plot true function

plot([0:.01:1], (constPoly+[0:.01:1].^expPoly(1)),'r');% true function

varIdStr = sprintf('%d', varId);
xaxis    = ['X', varIdStr];
yaxis    = ['Impact', varIdStr];
xlabel(xaxis,'FontSize',12)
ylabel(yaxis,'FontSize',12)
title(['Impact Function for ' xaxis],'FontSize',14)
legend('Estimated', 'True')

figure(2);
varId    = 2;%plot the intensity component for the second variable
bernsteinPolyPlot(beta,orderPoly,varId);%plot the estimated function
hold on;
plot([0:.01:1], [0:.01:1].^expPoly(2),'r');%true function

varIdStr = sprintf('%d', varId);
xaxis    = ['X', varIdStr];
yaxis    = ['Impact', varIdStr];
xlabel(xaxis,'FontSize',12)
ylabel(yaxis,'FontSize',12)
title(['Impact Function for ' xaxis],'FontSize',14)
legend('Estimated', 'True')
