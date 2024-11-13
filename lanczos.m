function rho = lanczos(A,v,lamCon,numeval)

% Lanczos method without restarting and reorth
% Only ritz values are returned.
if nargin==3
    numeval = 1;
end
v = v/norm(v);
iterMax = 4000;
T = zeros(iterMax,iterMax);
rho = inf*ones(numeval,iterMax);
for iter = 1:iterMax
    w = A(v);
    alpha = v' * w;
    if iter>1
        w = w-alpha*v-beta*vv;
    else
        w = w-alpha*v;
    end
    vv = v;
    [v,beta] = qr(w,0);
    T(iter,iter) = alpha;
    T(iter+1,iter) = beta;
    T(iter,iter+1) = beta;
    if iter>=numeval
        rho(:,iter) = eigs(T(1:iter,1:iter),numeval,'smallestabs');
    end
    if sum(rho(:,iter))<lamCon
        rho = rho(:,1:iter);
        break
    end
end


end
