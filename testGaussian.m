clear all
% close all

seedMax = 100;
histLC = cell(1,seedMax);

histKS = cell(1,seedMax);

nx = 300;
A = delsq(numgrid('L',2+nx));
n = size(A,1);
A = n*A;
lambda = eigs(A,1,'smallestabs');
numeval = 1;


dMax = 60;
iterMax = 1000;
tol = 1e-8;
toleval = lambda*(1+tol);
para.lamCon = toleval;
para.hist = 1;
para.iterMax = iterMax;
para.orth = 1;
para.tolra = 1e-6;
dMin = 30;

parfor iterSeed = 1:seedMax
    iterSeed
    rng(iterSeed);

    x0 = randn(n,1);
    [~,~,histLC{iterSeed}] = LC(@(x) A*x,x0,dMax,numeval,para);
    [~,~,histKS{iterSeed}] = lanczosTR(@(x) A*x,x0,dMax,numeval,dMin,para);

end




histLCs = cell(1,seedMax);
histKSs = cell(1,seedMax);
lambda = eigs(A,4,'smallestabs');
toleval = sum(lambda)*(1+tol);
numeval = 4;
para.lamCon = toleval;
x0 = randn(n,1);


parfor iterSeed = 1:seedMax
    iterSeed
    rng(iterSeed);

    x0 = randn(n,1);
    [~,~,histLCs{iterSeed}] = LC(@(x) A*x,x0,dMax,numeval,para);
    [~,~,histKSs{iterSeed}] = lanczosTR(@(x) A*x,x0,dMax,numeval,dMin,para);

end
















clear A
save('data_Gaussian')


jjMax = 5;
table1 = zeros(seedMax,jjMax);
table2 = zeros(seedMax,jjMax);
for ii = 1:seedMax
    for jj = 1:jjMax
        table1(ii,jj) = find((histKS{ii}.rho-lambda(1))/lambda(1)<10^(-3-jj),1);
        table2(ii,jj) = find((histLC{ii}.rho-lambda(1))/lambda(1)<10^(-3-jj),1);
    end
end

table = 1-table2./table1;
stat = genStat(table);
latex(sym(stat'))

jjMax = 5;
table1s = zeros(seedMax,jjMax);
table2s = zeros(seedMax,jjMax);
for ii = 1:seedMax
    for jj = 1:jjMax
        table1s(ii,jj) = find((sum(histKSs{ii}.rho)-sum(lambda))/sum(lambda)<10^(-3-jj),1);
        table2s(ii,jj) = find((sum(histLCs{ii}.rho)-sum(lambda))/sum(lambda)<10^(-3-jj),1);
    end
end

tables = 1-table2s./table1s;
stats = genStat(tables);
latex(sym(stats'))
clear A
save('data_Gaussian')


function stat = genStat(data)

[n,m] = size(data);
stat = zeros(4,m);
alpha = 0.05; 
for ii = 1:m
    outputs = data(:,ii);
    mean_value = mean(outputs);
    std_dev = std(outputs);  
    t_value = tinv(1 - alpha/2, n - 1); 
    sem = std_dev / sqrt(n);
    stat(1,ii) = mean_value;
    stat(2,ii) = std_dev;
    stat(3,ii) = mean_value - t_value * sem;
    stat(4,ii) = mean_value + t_value * sem;

    
end

end










