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





lanstep = 2; % minimal Lanczos step before restarting


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
            while norm(dz)>sqrt(j)*norm(w)*1e-15
                dz = V(:,1:j)'*w;
                w = w-V(:,1:j)*dz;
                z = z+dz;
            end
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
            errRitz(j+matvec-j0+1) = norm(errF);
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
                    [qs,~] = eigs(Ts(1:jj,1:jj),numeval,theta(1));
                    errF = beta*qs(jj,:);
                    vs = V(:,1:ii)*Qt(:,idxt(1:numeval));
                    zs = A(vs);
                    errapp(1,matvec-dMax+ii+1) = norm(zs-vs*(vs'*zs));
                    errapp(2,matvec-dMax+ii+1) = norm(errF);
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

        %mu=(1-sqrt(a/b))/(1+sqrt(a/b));
        %mSet(ii)=ii+ceil(log(2/para.tolra + 1) / (pi*ellipke(sqrt(1-mu^2))/(4*ellipke(mu))) - 1);
        mSet(ii)=ii+ceil(2*log(4/para.tolra)* log(4*b/a)/ (pi^2));
    end
    [m,ii] = min(mSet);
    m = min(m,dMax-lanstep)-ii;
    blksize(iter) = m+ii;


%     generate basis for rational Krylov subspace
    a = (theta(ii+1)-theta(numeval))/2;
    bb = (theta(ii+1)+theta(numeval))/2;
    b = max(Theta)-bb;


%     sq=sqrt(1-(a/b)^2);
%     K = ellipke(sq);
%     s = ceil(m/2);
%     poles=zeros(1,2*s);
%     for i=1:s
%        sn=ellipj((2*i-1)*K/(2*s),sq);
%        c=(a)*sqrt(-sn^2/(1-sn^2));
%        poles(2*i-1)=c;
%        poles(2*i)=-c;
%     end
%     poles = poles+bb;

    poles = zolopoles(a, b, m-1);
    poles = [poles+bb, inf];

    U = rat_krylov(T(1:dMax-1,1:dMax-1),b0,poles,'real');


    Q = orth([Q(:,1:ii),U]);
    j0 = size(Q,2)+1;
    b0 = Q'*b0;
    T = [Q'*T(1:dMax-1,1:dMax-1)*Q,b0;b0',0];
    T = (T+T')/2;
    T = blkdiag(T,zeros(dMax-j0,dMax-j0));
    V(:,1:j0-1) = V(:,1:dMax-1)*Q;
    V(:,j0) = v;

end



end


function xi = zolopoles(a,b, m)
    % xi is a vector containing the zolotarev poles in [-b,-a], [a,b] 
    kp = a/b;
    mu=(1-sqrt(kp))/(1+sqrt(kp));
    m_zol=ceil(m/2);
    
    alpha = acos(kp);
    K = mellipke(alpha);
    c = zeros(1,m_zol);
    
    for ii = 1:m_zol
       [sn,cn,~] = mellipj((2*ii-1)*K/(2*m_zol+1),alpha);   
       c(ii) = a * sn/cn; %we already compute the square root
    end
    
    xi = zeros(1, 2*m_zol);
    for i = 1:m_zol
        xi(2*i-1) = 1i*c(i);
        xi(2*i) = -1i*c(i);
    end
    
    end
    
    
    function [k,e] = mellipke(alpha,tol)
    %ELLIPKE Complete elliptic integral. Modified from Matlab's built-in code
    %for improved accuracy. 
    if nargin<1
      error(message('MATLAB:ellipke:NotEnoughInputs')); 
    end
    m = sin(alpha) * sin(alpha);
    m1 = cos(alpha) * cos(alpha);
    
    %classin = superiorfloat(m);
    classin='double';
    
    if nargin<2, tol = eps(classin); end
    if ~isreal(m),
        error(message('MATLAB:ellipke:ComplexInputs'))
    end
    if isempty(m), k = zeros(size(m),classin); e = k; return, end
    if any(m(:) < 0) || any(m(:) > 1), 
      error(message('MATLAB:ellipke:MOutOfRange'));
    end
    
    a0 = 1;
    b0 = cos(alpha);
    s0 = m;
    i1 = 0; mm = 1;
    while mm > tol
        a1 = (a0+b0)/2;
        b1 = sqrt(a0.*b0);
        c1 = (a0-b0)/2;
        i1 = i1 + 1;
        w1 = 2^i1*c1.^2;
        mm = max(w1(:));
        s0 = s0 + w1;
        a0 = a1;
        b0 = b1;
    end
    k = pi./(2*a1);
    e = k.*(1-s0/2);
    % im = find(m ==1);
    % if ~isempty(im)
    %     e(im) = ones(length(im),1);
    %    k(im) = inf;
    % end
    
    end
    
    
    function [sn,cn,dn] = mellipj(u,alpha,tol)
    %mELLIPJ Jacobi elliptic functions. MATLAB's built-in code, modified for improved accuracy. 
    %   [SN,CN,DN] = ELLIPJ(U,M) returns the values of the Jacobi elliptic 
    %   functions Sn, Cn and Dn, evaluated for corresponding elements of 
    %   argument U and parameter M.  U and M must be arrays of the same 
    %   size or either can be scalar.  As currently implemented, M is 
    %   limited to 0 <= M <= 1. 
    %
    %   [SN,CN,DN] = ELLIPJ(U,M,TOL) computes the elliptic functions to
    %   the accuracy TOL instead of the default TOL = EPS.  
    %
    
    if nargin<2
      error(message('MATLAB:ellipj:NotEnoughInputs')); 
    end
    
    u=real(u);
    alpha=real(alpha);
    
    m = sin(alpha) * sin(alpha);
    m1 = cos(alpha) * cos(alpha);
    
    %classin = superiorfloat(u,m);
    classin='double';
    if nargin<3, tol = eps(classin); end
    
    
    if ~isreal(u) || ~isreal(m)
        error(message('MATLAB:ellipj:ComplexInputs'))
    end
    
    if isscalar(m), m = m(ones(size(u))); end
    if isscalar(u), u = u(ones(size(m))); end
    if ~isequal(size(m),size(u)) 
      error(message('MATLAB:ellipj:InputSizeMismatch')); 
    end
    
    mmax = numel(u);
    
    %cn = zeros(size(u),classin);
    %cn = zeros(size(u));
    cn=u;
    sn = cn;dn = sn;
    
    m = m(:).';    % make a row vector
    u = u(:).';
    
    if any(m < 0) || any(m > 1), 
      error(message('MATLAB:ellipj:MOutOfRange'));
    end
    
    % pre-allocate space and augment if needed
    chunk = 10;
    %{
    a = zeros(chunk,mmax);
    c = a;
    b = a;
    a(1,:) = ones(1,mmax);
    c(1,:) = sin(alpha);
    b(1,:) = cos(alpha);
    %}
    c(1) = real(sin(alpha));
    b(1) = real(cos(alpha));
    a(1) = c(1);
          a = [a; zeros(chunk,mmax)];
          b = [b; zeros(chunk,mmax)];
          c = [c; zeros(chunk,mmax)];
          a(1)=1;
    
    n = zeros(1,mmax);
    i = 1;
    while any(abs(c(i,:)) > tol) && i<1000
        i = i + 1;
        if i > size(a,1)
          a = [a; zeros(chunk,mmax)];
          b = [b; zeros(chunk,mmax)];
          c = [c; zeros(chunk,mmax)];
        end
        a(i,:) = 0.5 * (a(i-1,:) + b(i-1,:));
        b(i,:) = sqrt(a(i-1,:) .* b(i-1,:));
        c(i,:) = 0.5 * (a(i-1,:) - b(i-1,:));
        in = find((abs(c(i,:)) <= tol) & (abs(c(i-1,:)) > tol));
        if ~isempty(in)
          [mi,ni] = size(in);
          n(in) = repmat((i-1), mi, ni);
        end
    end
    %phin = zeros(i,mmax,classin);
    phin(1)=u;
    phin = [phin;zeros(i-1,mmax)];
    
    
    phin(i,:) = (2 .^ n).*a(i,:).*u;
    while i > 1
        i = i - 1;
        in = find(n >= i);
        phin(i,:) = phin(i+1,:);
        if ~isempty(in)
    %      phin(i,in) = 0.5 * ...
    %      (asin(c(i+1,in).*sin(rem(phin(i+1,in),2*pi))./a(i+1,in)) + phin(i+1,in));
    %      phin(i,in) = 0.5 * ...
    %      (asin(c(i+1,in).*sin(rem(real(phin(i+1,in)),2*pi))./a(i+1,in)) + phin(i+1,in));
          phin(i,in) = 0.5 * ...
          (asin(c(i+1,in).*sin(real(phin(i+1,in)))./a(i+1,in)) + phin(i+1,in));
        end
    end
    %sn(:) = sin(rem(real(phin(1,:)),2*pi));
    %cn(:) = cos(rem(real(phin(1,:)),2*pi));
    sn(:) = sin(real(phin(1,:)));
    cn(:) = cos(real(phin(1,:)));
    dn(:) = sqrt(1 - m .* (sn(:).').^2);
    
    % special case m = 1 
    %m1 = find(m==1);
    %sn(m1) = tanh(u(m1));
    %cn(m1) = sech(u(m1));
    %dn(m1) = sech(u(m1));
    % special case m = 0
    %dn(m==0) = 1;
    
    
    
    end