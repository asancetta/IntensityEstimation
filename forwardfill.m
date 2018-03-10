function out = forwardfill(x)
[n,K]              = size(x);
out                = nan(n,K);
for k=1:K
   out(:,k) = forwardfillUtil(x(:,k));
    
end

end

function out  = forwardfillUtil(v)

idx = (~isnan(v)); %non nans
vr  = v(idx); %v non nan
cumIdx  = cumsum(idx);
nNan    = sum(cumIdx  == 0 ); 
out     = nan(size(v));
out(nNan+1:end) = vr(cumIdx(nNan+1:end)); %use cumsum to build index into vr


end

