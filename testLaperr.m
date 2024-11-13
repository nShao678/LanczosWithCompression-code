clear all
close all
rng(1);

nx = 300;
A = delsq(numgrid('L',2+nx));
n = size(A,1);
A = n*A;
lambda = eigs(A,4,'smallestabs');

x0 = randn(n,1);

dMax = 60;
iterMax = 5000;
tol = 1e-8;
para.hist = 1;
para.iterMax = iterMax;
para.orth = 1;
para.tolr = 1e-7;
para.F = 0;
para.tolra = 1e-4;
para.err = 1;


para.lamCon = lambda(1)*(1+tol);
numeval = 1;
[~,~,histLC1] = LC(@(x) A*x,x0,dMax,numeval,para);
rho1 = lanczos(@(x) A*x,x0,para.lamCon,numeval);

idx = find(rho1>para.lamCon);
figure
hold on
plot(histLC1.rho(1,idx)-lambda(1),'r-','LineWidth',2,'DisplayName','LC')
plot(rho1-lambda(1),'b--','LineWidth',2,'DisplayName','Lanczos')
plot(abs(histLC1.rho(1,idx)-rho1(idx)),'kx','LineWidth',1,'DisplayName','Difference')
hold off
axis([-inf,inf,-inf,inf])
set(gca,'yscale','log')
legend('FontSize',18,'Location','west')
set(gcf, 'Color', 'w');
export_fig('fig/exprho1.pdf')
export_fig('fig/exprho1.eps')

figure
hold on
plot(histLC1.errapp(1,:),'r-','LineWidth',2,'DisplayName','Norm of residual')
plot(histLC1.errapp(2,:),'b--','LineWidth',2,'DisplayName','Approximation')
plot(abs(histLC1.errapp(1,:)-histLC1.errapp(2,:)),'kx','LineWidth',2,'DisplayName','Difference')

hold off
axis([-inf,inf,-inf,inf])
set(gca,'yscale','log')
legend('FontSize',18,'Location','southwest')
set(gcf, 'Color', 'w');
export_fig('fig/experr1.pdf')
export_fig('fig/experr1.eps')

numeval = 4;
para.lamCon = sum(lambda(1:numeval))*(1+tol);

[~,~,histLC2] = LC(@(x) A*x,x0,dMax,numeval,para);
rho2 = lanczos(@(x) A*x,x0,para.lamCon,numeval);
idx = find(sum(rho2,1)>para.lamCon);
figure
hold on
plot(sum(histLC2.rho,1)-sum(lambda(1:numeval)),'r-','LineWidth',2,'DisplayName','LC')
plot(sum(rho2,1)-sum(lambda(1:numeval)),'b--','LineWidth',2,'DisplayName','Lanczos')
plot(abs(sum(rho2(:,idx),1)-sum(histLC2.rho(:,idx),1)),'kx','LineWidth',2,'DisplayName','Difference')

hold off
axis([-inf,inf,-inf,inf])
set(gca,'yscale','log')
legend('FontSize',18,'Location','west')
set(gcf, 'Color', 'w');
export_fig('fig/exprho2.pdf')
export_fig('fig/exprho2.eps')

figure
hold on
plot(histLC2.errapp(1,:),'r-','LineWidth',2,'DisplayName','Norm of residual')
plot(histLC2.errapp(2,:),'b--','LineWidth',2,'DisplayName','Approximation')
plot(abs(histLC2.errapp(1,:)-histLC2.errapp(2,:)),'kx','LineWidth',2,'DisplayName','Difference')


hold off
axis([-inf,inf,-inf,inf])
set(gca,'yscale','log')
legend('FontSize',18,'Location','southwest')
set(gcf, 'Color', 'w');
export_fig('fig/experr2.pdf')
export_fig('fig/experr2.eps')

clear A
save('data_testLaperr')







