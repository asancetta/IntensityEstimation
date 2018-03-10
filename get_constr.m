function [A, c] = get_constr(indVar,k,a)
% indVar index that identifies all the columns that belong to the same
% variable. For example the data is such that 
% data = [x(:,1),x(:,1).^2, x(:,1).^3, x(:,2), x(:,2)], then indVar = [1,1,1,2,2].
% k is the variable we are interested in, e.g. k=1, first variable
% a is the restriction, for example, a = [1,-1] for monotonic increasing, 
% a = -[1,-1] for monotonic decreasing, a = [-1,2,-1] for convexity
% the restriction is data*A < zero. 
nA                                = length(a)-1;
K1                                = length(indVar);
ind1                              = find(indVar == k);
k1                                = length(ind1);
A                                 = zeros(k1-nA,K1);
for i = 1:(k1-nA)
    
    A(i,ind1(i:i+nA))             = a;

end
c                                 = zeros(k1-nA,1);


end
