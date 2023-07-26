function [W,V]=GCFG(X,Z_ini,lambda1,max_iter,miu)
%%
max_miu = 1e8;     
tol  = 1e-6;
tol2 = 1e-2;
rho = 3;
[m,n] = size(X);
C1 = zeros(n,n);
V_ini=eye(n);

%%
    % ----------Initilization for W, V--------% 

for iter = 1:max_iter
    if iter == 1
        W = Z_ini;
        V = V_ini;
    end
    W_old = W;
    V_old = V;
    
    % ---------- Optimization for W  -------- %
    Ctg = X'*X*W*V'*V+lambda1*W;
    W = Ctg\(W*X'*X*V);
    
    % ---------- Optimization for V  -------- %
    G1 = W'*X'*X*W+C1;
     V = (X'*X*W)/(G1);    

    % ------ C1 C2 miu ---------- %
    L1 = V'*V-eye(size(V,2));
    C1 = C1+miu*L1;
    
    LL1 = norm(W-W_old,'fro');
    LL2 = norm(V-V_old,'fro');
    SLSL = max(max(LL1,LL2))/norm(X,'fro');
    if miu*SLSL < tol2
        miu = min(rho*miu,max_miu);
    end
    stopC = norm(L1,'fro');
    if stopC < tol
        iter
        break;
    end
    obj(iter) = stopC;
end

