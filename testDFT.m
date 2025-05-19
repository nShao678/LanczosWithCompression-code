clear all
close all
rng(1);



name = {'H2O','GaAsH6','SiO2','Si5H12','Ga10As10H30','Si34H36','Si41Ge41H72','Si87H76','Ge99H100'};
len = length(name);
histLC = cell(1,len);
histKS = cell(1,len);
numevalSet = [4,7,8,16,55,86,200,212,248];
matvecSet = [500,400,600,600,1500,1500,3500,3500,3500];
matvec = zeros(2,len);
rhoerr = zeros(2,len);
result = zeros(4,len);
len = 9;
for ii = 1:len
    load(['data/',name{ii},'.mat']);
    A = Problem.A;
    n = size(A,1);


    numeval = numevalSet(ii);
    dMax = max(4*numeval,80);
    lambda = eigs(A,numeval,'smallestreal','MaxIterations',100,'SubspaceDimension',dMax);
    x0 = randn(n,1);


    para.hist = 1;
    para.matvecMax = matvecSet(ii);
    para.tolra = 1e-8;
    para.F = 0;
    para.orth = 1;
    dMin = ceil(dMax/2);
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

name = {'H2O','GaAsH6','SiO2','Si5H12','Ga10As10H30','Si34H36','Si41Ge41H72','Si87H76','Ge99H100'};
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
    export_fig(['fig/',name{ii},'.pdf'])
    export_fig(['fig/',name{ii},'.eps'])
end


