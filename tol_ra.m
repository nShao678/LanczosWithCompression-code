clear all
close all
rng(1);
iiMax = 5;
histLC = cell(1,iiMax);
histKS = cell(1,iiMax);
histLan = cell(1,iiMax);
lambdaSet = zeros(1,iiMax);
errSet = cell(1,iiMax);

numeval = 1;
for ii = 1:iiMax
    
    tic
    nx = 300;
    A = delsq(numgrid('L',2+nx));
    n = size(A,1);
    A = n*A;
    lambda = sum(eigs(A,numeval,'smallestabs'));
    x0 = randn(n,1);

    dMax = 60;
    tol = 0;
    toleval = lambda*(1+tol);
    para.lamCon = toleval;
    para.hist = 1;
    para.matvecMax = 1500;
    para.orth = 1;
    para.tolra = 10^(-ii);

    [~,~,histLC{ii}] = LC(@(x) A*x,x0,dMax,numeval,para);
    errSet{ii} = (sum(histLC{ii}.rho,1)-lambda)./lambda;
    toc
    ii
end
save('data_tol_ra_s')


numeval = 4;
for ii = 1:iiMax
    
    tic
    nx = 300;
    A = delsq(numgrid('L',2+nx));
    n = size(A,1);
    A = n*A;
    lambda = sum(eigs(A,numeval,'smallestabs'));
    x0 = randn(n,1);

    dMax = 60;
    tol = 0;
    toleval = lambda*(1+tol);
    para.lamCon = toleval;
    para.hist = 1;
    para.matvecMax = 1500;
    para.orth = 1;
    para.tolra = 10^(-ii);

    [~,~,histLC{ii}] = LC(@(x) A*x,x0,dMax,numeval,para);
    errSet{ii} = (sum(histLC{ii}.rho,1)-lambda)./lambda;
    toc
    ii
end
save('data_tol_ra_m')


%%
load("data_tol_ra_s.mat")
figure

hold on
for ii = 1:iiMax
plot(errSet{ii},'linewidth',2,'DisplayName',['$\mathsf{tol}_{\ \mathrm{ra}}=',num2str(10^(-ii),'$%0.0e')])
end
hold off
ylabel('Relative error of Ritz values')
legend('FontSize',18,'Location','southwest','Box','off','Interpreter','latex')
set(gcf, 'Color', 'w');
set(gca,'yscale','log');
export_fig('fig/LapTolraSingle.pdf')
export_fig('fig/LapTolraSingle.eps')


load("data_tol_ra_m.mat")
figure

hold on
for ii = 1:iiMax
plot(errSet{ii},'linewidth',2,'DisplayName',['$\mathsf{tol}_{\ \mathrm{ra}}=',num2str(10^(-ii),'$%0.0e')])
end
hold off
ylabel('Relative error of Ritz values')
legend('FontSize',18,'Location','southwest','Box','off','Interpreter','latex')
set(gcf, 'Color', 'w');
set(gca,'yscale','log');
export_fig('fig/LapTolraSeveral.pdf')
export_fig('fig/LapTolraSeveral.eps')

