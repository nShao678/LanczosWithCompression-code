function [u,matvec,hist] = LC(A,v,dMax,numeval,para)
% dMax: The maximum allowed dimension of Krylov space.
% numeval: The number of desired eigenvalues
% lamCon: The stopping criteria is ritz value smaller than lamCon.
% full reorthogonalization is applied


if nargin==4
    para = [];
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
if ~isfield(para,'orth')
    para.orth = 1;      % reorth will fill-in (recommand)
end
if ~isfield(para,'hist')
    para.hist = 0;      % convergence history (expensive)
end
if ~isfield(para,'err')
    para.err = 0;       % history of residual (expensive)
end
if ~isfield(para,'tolra')
    para.tolra = 1e-12; % rational approximation
end
if ~isfield(para,'checkRitz')
    para.checkRitz = 0; % check error in Ritz (expensive)
end





if para.err == 1
    para.hist = 1;
    errapp = zeros(2,para.matvecMax+dMax);
end
if para.hist == 1
    para.tolr = -inf;
else
    para.lamCon = -inf;
end




lanstep = 5; % minimal Lanczos step before restarting


n = size(v,1);
V = zeros(n, dMax);
T = zeros(dMax, dMax);

v = orth(v);
V(:, 1) = v;
nA = 0;
j0 = 1;

rho = inf*ones(numeval,para.matvecMax+dMax);

matvec = 0;
blksize = zeros(1,para.matvecMax+dMax);
tau = cell(1,para.matvecMax+dMax);

if para.err==1 || para.tolr ~= -inf
    alphaSet = zeros(1,para.matvecMax+dMax);
    betaSet = zeros(1,para.matvecMax+dMax);
    Ts = [];
end
errF = 0;
errRitz = zeros(1,para.matvecMax);


for iter = 1:para.matvecMax
    for j = j0:dMax-1
        w = A(v);
        % sacrifices tridiagonal structure of $T$ to improve stability
        if para.orth==1
            z = V(:,1:j)'*w;
            w = w-V(:,1:j)*z;
            dz = V(:,1:j)'*w;
            w = w-V(:,1:j)*dz;
            z = z+dz;
            dz = V(:,1:j)'*w;
            w = w-V(:,1:j)*dz;
            z = z+dz;
            T(1:j,j) = z;
            T(j,1:j-1) = z(1:j-1)';
            alpha = T(j,j);
            [v,beta] = qr(w,0);
            V(:,j+1) = v;
        else
            alpha = v' * w;
            T(j,j) = alpha;
            w = w-V(:,1:j)*(V(:,1:j)'*w);
            w = w-V(:,1:j)*(V(:,1:j)'*w);
            [v,beta] = qr(w,0);
            V(:,j+1) = v;
            T(j,j+1) = beta;
            T(j+1,j) = beta;
        end
        if para.checkRitz == 1 
            AV = A(V(:,1:j));
            errF = V(:,1:j)'*AV-T(1:j,1:j);
            errRitz(j+matvec-j0+1) = norm(errF,'fro');
        end
        if para.err==1 || para.tolr ~= -inf
            alphaSet(matvec+j-j0+1) = alpha;
            betaSet(matvec+j-j0+1) = beta;
        end
    end

    [Q,Theta] = eig(T(1:dMax-1,1:dMax-1),'vector');
    [theta,idx] = mink(Theta,dMax-1-lanstep);
    matvec0 = matvec;
    matvec = matvec+dMax-j0;
    if para.err == 1 || para.tolr ~= -inf
        nA = max([abs(Theta);nA]);
        Ts = diag(alphaSet(1:matvec))+diag(betaSet(1:matvec-1),1)+diag(betaSet(1:matvec-1),-1);
    end
    if para.hist == 1
        % To plot the convergence history
        for ii = j0:dMax-1
            if ii>=numeval
                [Qt,Dt] = eig(T(1:ii,1:ii),'vector');
                [rho(:,matvec-dMax+ii+1),idxt] = mink(Dt,numeval);
                VV = orth(V(:,1:ii)*Qt(:,idxt));
                rho(:,matvec-dMax+ii+1) = sort(eig(VV'*A(VV),'vector'));
                if para.err == 1
                    jj = matvec0+ii-j0+1;
                    [qs,~] = eigs(Ts(1:matvec,1:matvec),numeval,theta(1));
                    errF = beta*qs(jj,:);
                    vs = V(:,1:ii)*Qt(:,idxt(1:numeval));
                    zs = A(vs);
                    errapp(1,matvec-dMax+ii+1) = norm(zs-vs*(vs'*zs),'fro');
                    errapp(2,matvec-dMax+ii+1) = norm(errF,'fro');
                end
            end
        end
    elseif para.err == 1 || para.tolr ~= -inf
        [qs,~] = eigs(Ts(1:matvec,1:matvec),numeval,theta(1));
        errF = beta*qs(matvec,:);
    end

    if sum(rho(:,matvec))<para.lamCon || matvec>=para.matvecMax || norm(errF)/nA<para.tolr
        u = V(:,1:dMax-1)*Q(:,idx(1:numeval));
        if para.hist == 1
            hist.rho = rho(:,1:matvec);
            if para.err==1
                hist.errapp = errapp(:,1:matvec);
            end
        end
        if para.checkRitz == 1
            hist.errRitz = errRitz(1:matvec);
        end
        hist.blksize = blksize(1:iter-1);
        hist.tau = tau(1:iter-1);
        break
    end

    b0 = zeros(dMax-1,1);
    b0(dMax-1,:) = beta';
    Q = Q(:,idx);


    % find the position of rational function
    mSet = inf*ones(1,dMax-2-lanstep);
    for ii = numeval:length(mSet)
        a = (theta(ii+1)-theta(numeval))/2;
        bb = (theta(ii+1)+theta(numeval))/2;
        b = max(Theta)-bb;

        mu=(1-sqrt(a/b))/(1+sqrt(a/b));
        mSet(ii)=ii+ceil(log(2/para.tolra + 1) / (pi*ellipke(sqrt(1-mu^2))/(4*ellipke(mu)))-1);
    end
    [m,ii] = min(mSet);
    m = min(m,dMax-lanstep)-ii;
    blksize(iter) = m+ii;


    % generate basis for rational Krylov subspace
    a = (theta(ii+1)-theta(numeval))/2;
    bb = (theta(ii+1)+theta(numeval))/2;
    tau{iter}.tau = bb;
    tau{iter}.theta = theta;
    b = max(Theta)-bb;
    sq=sqrt(1-(a/b)^2);
    K = ellipke(sq);
    s=ceil(m/2);
    poles=zeros(1,2*s);
    for i=1:s
        sn=ellipj((2*i-1)*K/(2*s),sq);
        c=(a)*sqrt(-sn^2/(1-sn^2));
        poles(2*i-1)=c;
        poles(2*i)=-c;
    end
    poles = poles+bb;
    poles = [poles,inf];
    U = rat_krylov(T(1:dMax-1,1:dMax-1),b0,poles,'real');
    Q = Q(:,1:ii);
    U = orth(U-Q*(Q'*U));
    U = orth(U-Q*(Q'*U));
    Q = [Q,U];
    j0 = size(Q,2)+1;
    b0 = Q'*b0;
    T = [Q'*T(1:dMax-1,1:dMax-1)*Q,b0;b0',0];
    T = (T+T')/2;
    T = blkdiag(T,zeros(dMax-j0,dMax-j0));
    V(:,1:j0-1) = V(:,1:dMax-1)*Q;
    V(:,j0) = v;

end



end
