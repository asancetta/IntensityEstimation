function out = weibullMeasure(time, isJump,beta,flagDensity, g)
% time and isJump have the same size, isJump is 1 if time is a jump time,
% zero otherwise, beta is the shape parameter in Weibull hazard,
% flagDensity =0 if we want to compute the integral 


if nargin  < 4
    flagDensity  = 0;
end

if nargin  < 5
   g  = 1;
end


if flagDensity == 0 
    %compute the integral of the Weibull hazard over the observed times;
    %these times are not necessarily jump times
    jumpTimes          = makeItLeftContinuous(time,time,isJump);
    jumpTimes(1)       = 0; %assume there is a jump at time 0
    timeLag            = [nan(1,1); time(1:end-1)];
    out                = (time-jumpTimes).^beta  - (timeLag-jumpTimes).^beta;
    out(out<0)         = 0;
    out                = out.*g;
    out(isnan(out))    = 0;
else
    %compute the Weibull hazard at jump times only
    out                = get_diff1(time(isJump));
    
    
    out                = (beta*out.^(beta-1));
    out(1)             = 1;
    out                = out.*g;

end



end


function out =  get_diff1(x)
[n,K]              = size(x);
out                = nan(n,K);
out(2:end,:)       = x(2:end,:)-x(1:end-1,:);

end

