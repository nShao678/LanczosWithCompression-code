clear all
close all
rng(1);



name = {'H2O','GaAsH6','SiO2','Si5H12','Ga10As10H30','Si34H36'};
len = length(name);
histLC = cell(1,len);
histKS = cell(1,len);
numevalSet = [4,7,8,16,55,86];
matvecSet = [500,500,600,600,1500,1500];
matvec = zeros(2,len);
rhoerr = zeros(2,len);
result = zeros(4,len);

for ii = 1:len
    load(['data/',name{ii},'.mat']);
    A = Problem.A;
    n = size(A,1);


    numeval = numevalSet(ii);
    dMax = max(4*numeval,80);
    [u,lambda] = eigs(A,numeval,'smallestreal','MaxIterations',100,'SubspaceDimension',dMax);
    lambda = diag(lambda);
    para.eval = lambda;
    x0 = randn(n,1);


    para.hist = 1;
    para.matvecMax = matvecSet(ii);
    para.tolra = 1e-7;
    para.F = 0;
    para.orth = 1;
    dMin = ceil(dMax/2);
    para.u0 = u;
    [uLC,matvec(1,ii),histLC{ii}] = LC(@(x) A*x,x0,dMax,numeval,para);
    histLC{ii}.rho(1:numeval,:) = histLC{ii}.rho(1:numeval,:)-lambda;

    [uKS,matvec(2,ii),histKS{ii}] = lanczosTR(@(x) A*x,x0,dMax,numeval,dMin,para);
    histKS{ii}.rho(1:numeval,:) = histKS{ii}.rho(1:numeval,:)-lambda; 
    figure
    hold on

    plot(sum(histLC{ii}.rho,1),'b-','linewidth',2,'DisplayName','LC')
    plot(sum(histKS{ii}.rho,1),'r--','linewidth',2,'DisplayName','KS')
    hold off
    legend('FontSize',18)
    set(gca,'yscale','log')
    clear A
    save('data_DFT')
    pause(1)
end
%%
name = {'H2O','GaAsH6','SiO2','Si5H12','Ga10As10H30','Si34H36'};
len = length(name);
for ii = 1:len
figure
    hold on

    plot(sum(histLC{ii}.rho,1),'b-','linewidth',2,'DisplayName','LC')
    plot(sum(histKS{ii}.rho,1),'r--','linewidth',2,'DisplayName','KS')
    hold off
    % legend('FontSize',18)
    set(gca,'yscale','log')
    set(gcf, 'Color', 'w');
    axis([-inf,inf,1e-14,inf])
    export_fig(['fig/',name{ii},'.pdf'])
    export_fig(['fig/',name{ii},'.eps'])
end


