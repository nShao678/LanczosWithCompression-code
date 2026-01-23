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