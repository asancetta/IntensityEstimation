function out = bernsteinPolyPlot(beta,orderPoly,varId)



nBeta   = length(beta);
K       = nBeta/(orderPoly+1);
indVar  =  kron((1:K),ones(1,orderPoly+1));
indActive   = indVar == varId;
betaActive  = beta(indActive);
    
x              = [0:.01:1]';
basis          = get_bernsteinBasisMultivariate(orderPoly, x);
out            = basis*betaActive;
plot(x, out)

end