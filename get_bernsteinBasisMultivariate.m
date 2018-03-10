function basis = get_bernsteinBasisMultivariate(orderPoly,X)


basis = arrayfun(@(i) bernsteinBasisUnivariate(orderPoly,X(:,i)),...
                 1:size(X,2), 'UniformOutput', false);
basis = cell2mat(basis);             



end


function poly  = bernsteinBasisUnivariate(m,x)
n              = size(x,1);
poly           = zeros(n,m+1);

for j=0:m

    poly(:,j+1)       = x.^(j).*(1-x).^(m-j)*nchoosek(m,j);    
    
end

end


