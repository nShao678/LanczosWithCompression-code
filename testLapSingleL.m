clear all
close all
rng(1);
nx = 300;
A = delsq(numgrid('L',2+nx));
n = size(A,1);
A = n*A;
lambda = eigs(A,1,'smallestabs');

x0 = randn(n,1);

dMax = 60;
iterMax = 10000;
tol = 1e-8;
toleval = lambda*(1+tol);
numeval = 1;
para.lamCon = toleval;
para.hist = 1;
para.matvecMax = iterMax;
para.orth = 1;
para.tolr = 1e-7;
para.F = 0;
para.tolra = 1e-4;

[~,~,histLC] = LC(@(x) A*x,x0,dMax,numeval,para);
dSet = [1,2,4,8,16,25,30,35,40];
iiMax = length(dSet);
histKS = cell(1,iiMax);
for ii = 1:iiMax
    dMin = dSet(ii)
    [~,~,histKS{ii}] = lanczosTR(@(x) A*x,x0,dMax,numeval,dMin,para);
end


figure
hold on
for ii = 1:iiMax
plot((histKS{ii}.rho-lambda)/lambda,'--','DisplayName',['KS-',num2str(dSet(ii))],'LineWidth',2)
end
plot((histLC.rho-lambda)/lambda,'b--','DisplayName','LC','LineWidth',2)
hold off
set(gca,'yscale','log')
legend('FontSize',18)
set(gcf, 'Color', 'w');
axis([-inf,inf,tol,inf])
export_fig('fig/expLapSingleL.pdf')
export_fig('fig/expLapSingleL.eps')

jjMax = 5;
table = zeros(iiMax+2,jjMax);
for ii = 1:iiMax
    for jj = 1:jjMax
        table(ii,jj) = find((histKS{ii}.rho-lambda)/lambda<10^(-3-jj),1);
    end
end
for jj = 1:jjMax
    table(iiMax+1,jj) = find((histLC.rho-lambda)/lambda<10^(-3-jj),1);
end
table(iiMax+2,:) = 1-table(iiMax+1,:)./min(table(1:iiMax,:));
latex(sym(table))


save('data_testLapSingleL')




