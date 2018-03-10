function out  = makeItLeftContinuous(X,time,isJump)
%time is the time at which X is observed, time(isJump) is the time when X
%changes; this assumes X is right continuous. The function makes X left
%continuous, i.e. if T_i's are jump times, shifts X rows so that
%X(t)=X(T_{i}) if t in (T_{i},T_{i+1}]

[n,K]                    = size(X);
epsTime                  = min(diff(time))/10;
timePlus                 = time;
timePlus(isJump)         = time(isJump)+epsTime;
timeAll                  = union(time,timePlus);
[timeMask, indTimePlus]  = intersect(timeAll,timePlus);
[timeMask, indTime]      = intersect(timeAll,time);

Xmask                    = nan(size(timeAll,1),K);
Xmask(indTimePlus,:)     = X;
Xmask                    = forwardfill(Xmask);
out                      = Xmask(indTime,:);

end


