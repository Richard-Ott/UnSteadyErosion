function [K,f] = ebc(K,f,bound,val)

%  function [K,f] = ebc(K,f,bound,val)
%
%    Enforces essential boundary conditions on nodes in vector bound
%    If val is a vector of length equal to the length of bound then 
%    each entry in this vector will be enforced on the associated node 
%    of vector bound. That is val(1) -> bound(1) ... val(n) -> bound(n)
%
%    In cases where the length of val differs from the length of bound
%    the first entry in val will be enforced on all nodes in bound.
%    That is val(1) -> bound(1..n)
%    As a special case, if a constant value is enforced on the entire
%    boundary val may be a scalar.
%
%    K and f are the global stiffness matrix and load vector, respectively.
%
%    See also nbc
%

%Loop nodes in bound
for i=1:length(bound),
    
    if size(bound) == size(val), lval = val(i);
    else, lval = val(1); end;
    
    %Modify load vector
    f = f - lval*K(:,bound(i));
    f(bound(i)) = lval;
    
    %Modify stiffness matrix
    K(:,bound(i)) = zeros(size(K(:,bound(i))));
    K(bound(i),:) = zeros(size(K(bound(i),:)));
    K(bound(i),bound(i)) = 1;
    
end;