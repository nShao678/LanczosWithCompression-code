clear all
close all
warning off
rng(1);



    load('data\Ga10As10H30.mat');
    A = Problem.A;
    n = size(A,1);
    numeval = 55;
    dMax = 4*numeval;
    lambda = eigs(A,numeval,'smallestreal','MaxIterations',100,'SubspaceDimension',dMax);
    


    x0 = randn(n,1);
    para.tolra = 1e-12;
    para.orth = 1;
    para.tolr = -inf;
    
    para.matvecMax = 2000;
    para.checkRitz = 1;
    para.hist = 0;
    para.orth = 0;
    para.err = 0;
    [~,~,hist0] = LC(@(x) A*x,x0,dMax,numeval,para);

    para.orth = 1;
    [~,~,hist1] = LC(@(x) A*x,x0,dMax,numeval,para);
    clear A
    save('data_DFTerr')

    figure
    hold on
    plot(hist1.errRitz,'b-','LineWidth',2,'DisplayName','With fill-in')
    plot(hist0.errRitz,'r--','LineWidth',2,'DisplayName','Without fill-in')
    hold off
    legend('FontSize',18,'Location','east')
    set(gca,'yscale','log')
    set(gcf, 'Color', 'w');
    axis([0,1600,-inf,inf])
    export_fig(['fig/reorth.pdf'])
    export_fig(['fig/reorth.eps'])

    


