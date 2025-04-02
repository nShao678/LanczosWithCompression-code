clear all
close all
rng(1);
iiMax = 10;
histLC = cell(1,iiMax);
histKS = cell(1,iiMax);
histLan = cell(1,iiMax);
lambdaSet = zeros(1,iiMax);
for ii = 1:iiMax
    
    tic
    nx = ii*100;
    A = delsq(numgrid('L',2+nx));
    n = size(A,1);
    A = n*A;
    lambda = eigs(A,4,'smallestabs');
    lambdaSet(ii) = sum(lambda);
    x0 = randn(n,1);

    dMax = 60;
    tol = 1e-8;
    toleval = lambdaSet(ii)*(1+tol);
    numeval = 4;
    para.lamCon = toleval;
    para.hist = 1;
    para.matvecMax = 6000;
    para.orth = 1;
    para.tolra = 1e-6;

    [~,~,histLC{ii}] = LC(@(x) A*x,x0,dMax,numeval,para);

    dMin = ceil(dMax/2);
    [~,~,histKS{ii}] = lanczosTR(@(x) A*x,x0,dMax,numeval,dMin,para);
    histLan{ii} = lanczos(@(x) A*x,x0,para.lamCon);
    toc
    ii
end


jjMax = 5;
table1 = zeros(iiMax,jjMax);
table2 = table1; 
table3 = table2;
for ii = 1:iiMax
    for jj = 1:jjMax
        table1(ii,jj) = find((sum(histLC{ii}.rho,1)-lambdaSet(ii))/lambdaSet(ii)<10^(-3-jj),1);
        table2(ii,jj) = find((sum(histKS{ii}.rho,1)-lambdaSet(ii))/lambdaSet(ii)<10^(-3-jj),1);
    end
end


locx = (1:iiMax).^2*7500;
figure
hold on
for jj = 1:jjMax
    plot(locx,100*(1-table1(:,jj)./table2(:,jj)),'x-','linewidth',2,'DisplayName',['relerr=',num2str(10^(-3-jj),'%0.0e')])
end
hold off
legend('FontSize',18,'Location','southeast','Box','off')
xlabel('Matrix size')
ytickformat('percentage');
ylabel('Improvement of LC over KS')
set(gcf, 'Color', 'w');
export_fig('fig/LapServeralN.pdf')
export_fig('fig/LapServeralN.eps')

latex(sym([locx',table1]))
latex(sym([locx',table2]))
latex(sym([locx',1-table1./table2]))


clear A
save('data_testLapSeveralN')




