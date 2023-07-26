function [ Z ] = GGM( X, lambda, dig )

if (dig)
    D = pinv(X' * X + lambda*eye(size(X,2)));
    Z = -lambda*D / spdiags(diag(D), 0, size(X,2), size(X,2));     
    Z(logical(eye(size(Z)))) = 0;                          
else
    Z = -lambda*(X' * X + lambda*eye(size(X,2))) \ X'*X;
end

