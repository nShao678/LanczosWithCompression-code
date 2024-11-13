function [u,matvec,hist] = lanczosTR(A,v,dMax,numeval,dMin,para)

%dMax: The maximum allowed dimension of Krylov space 
%dMin: The number of Ritz vectors kept after restarting 
%lamCon: The stopping criteria is ritz value smaller than tol
%full reorthogonalization is applied

if nargin==5
    para = []
end
if ~isfield(para,'matvecMax')
    para.matvecMax = 2000; % restarting cycle
end
if ~isfield(para,'lamCon')
    para.lamCon = -inf; % SC: sum(\mu_{i})<lamCon
end
if ~isfield(para,'tolr')
    para.tolr = 1e-7;   % SC: norm(R,fro)/nA<tolr
end
if ~isfield(para,'hist')
    para.hist = 0;      % convergence history (expensive)
end
if para.hist==1
    para.tolr = -inf;
end

[n,~] = size(v);
V = zeros(n, dMax);
T = zeros(dMax, dMax);
v = orth(v);
V(:, 1) = v;
nA = 0;
j0 = 1;
if para.hist == 1
    rho = zeros(numeval,para.matvecMax+dMax);
    err = zeros(numeval,para.matvecMax+dMax);
end
matvec = 0;
hist = [];

for iter = 1:para.matvecMax
    for j = j0:dMax-1
        w = A(v);
        alpha = v' * w;
        T(j,j) = alpha;
        w = w-V(:,1:j)*(V(:,1:j)'*w);
        w = w-V(:,1:j)*(V(:,1:j)'*w);
        [v,beta] = qr(w,0);
        V(:, j+1) = v;
        T(j, j+1) = beta;
        T(j+1,j) = beta;
    end

    [Q,Theta] = eig(T(1:dMax-1,1:dMax-1),'vector');
    nA = max([abs(Theta);nA]);
    [~,idx] = mink(Theta,dMin);
    matvec = matvec+dMax-j0;

    if para.hist==1
        for ii = j0:dMax-1
            if ii<numeval
                rho(:,matvec-dMax+ii+1) = inf;
                err(:,matvec-dMax+ii+1) = inf;
            else
                [Qt,Dt] = eig(T(1:ii,1:ii),'vector');
                [rho(:,matvec-dMax+ii+1),idxt] = mink(Dt,numeval);
                for jj = 1:numeval
                    err(jj,matvec-dMax+ii+1) = norm(T(ii+1,ii)*Qt(ii,idxt(jj)),'fro')/nA;
                end
            end
        end
    end


    if sum(mink(Theta,numeval))<para.lamCon || matvec>=para.matvecMax || norm(T(dMax, dMax-1)*Q(dMax-1,idx(1:numeval)),'fro')/nA<para.tolr
        u = V(:,1:dMax-1)*Q(:,idx(1));
        if para.hist == 1
            hist.rho = rho(:,1:matvec);
            hist.err = err(:,1:matvec);
        end
        break
    end

    V(:,1:dMin) = V(:,1:dMax-1)*Q(:,idx);
    V(:,dMin+1) = v;
    T = [diag(Theta(idx)),Q(dMax-1,idx)'*beta;beta*Q(dMax-1,idx),0];
    j0 = dMin+1;
end


end
