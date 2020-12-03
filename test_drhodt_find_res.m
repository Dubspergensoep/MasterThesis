%SOLVE-DRHODT-FIND_RES test
% system/hamiltonian parameters
gamma=0.005;
Kappa=0.005;
g=0.5;
J=0.126;
omega=1;
eps=omega-J;
mu=0.75;
% sym parameters
fill=0;
spin=[0;1];
N_max=5;
s_WF=2*(N_max+1);       %size wave function
sig = GetAllOperatorsSparse(1,N_max);
% Laser parameters
fmax=0.1;          %-> dit moet nog bepaald worden
rate=1e-6;
twait=1000;
fmin=0.05;
resstep=10;
extra_steps=100;
method.str='der_abs';
method.threshold=30;
method.abs_threshold=0.1;

% reduction parameters
dT=100;
% Construct rho_0
Psi=GetFockWf(sig,1,fill,spin);
Rho_0=Psi*Psi';
% Construct H1
H1=H_MF_simple1(sig,omega,mu,eps,g);
%% Solve
[rho,timecell,fcell]=Solve_dhrodt_find_res(H1,sig,J,gamma,Kappa,Rho_0,dT,rate,fmin,fmax,twait,resstep,extra_steps,method);
%% expectation values
v_ea_inc=expectation_val(length(timecell{1}),rho{1},sig{1,2});
v_ea_wait=expectation_val(length(timecell{2}),rho{2},sig{1,2});
v_ea_dec=expectation_val(length(timecell{3}),rho{3},sig{1,2});
v_ea=[v_ea_inc v_ea_wait v_ea_dec];
time=[timecell{1} timecell{2} timecell{3}];
%% plots
plottitle=['JCH bistability: J=' num2str(J) ' \gamma=\kappa=' num2str(gamma) '  \mu=' num2str(mu) '  \epsilon=\omega - J'];
saveplot.saveplot=false;
vf=[fcell{1} fcell{2} fcell{3}];
plot_a2f2(time,v_ea,vf,plottitle,saveplot)
indices=[length(v_ea_inc) (length(v_ea_inc)+length(v_ea_wait))];
plot_a2f2color(v_ea,vf,plottitle,saveplot,indices)
plot_a2f2color2(v_ea,time,plottitle,saveplot,indices)
%% A NON RESONANT EXAMPLE
mu=0.3;
g=0.7;
fmax=0.5
% Construct H1
H1=H_MF_simple1(sig,omega,mu,eps,g);
%% Solve
[rho,timecell,fcell]=Solve_dhrodt_find_res(H1,sig,J,gamma,Kappa,Rho_0,dT,rate,fmin,fmax,twait,resstep,extra_steps,method);
%% expectation values
v_ea_inc=expectation_val(length(timecell{1}),rho{1},sig{1,2});
v_ea_wait=expectation_val(length(timecell{2}),rho{2},sig{1,2});
v_ea_dec=expectation_val(length(timecell{3}),rho{3},sig{1,2});
v_ea=[v_ea_inc v_ea_wait v_ea_dec];
time=[timecell{1} timecell{2} timecell{3}];
%% plots
plottitle=['JCH bistability: J=' num2str(J) ' \gamma=\kappa=' num2str(gamma) '  \mu=' num2str(mu) '  \epsilon=\omega - J'];
saveplot.saveplot=false;
vf=[fcell{1} fcell{2} fcell{3}];
plot_a2f2(time,v_ea,vf,plottitle,saveplot)
indices=[length(v_ea_inc) (length(v_ea_inc)+length(v_ea_wait))];
plot_a2f2color(v_ea,vf,plottitle,saveplot,indices)
plot_a2f2color2(v_ea,time,plottitle,saveplot,indices)

%% Check new find_bistab function
mu=0.75;
g=0.5;
fmax=0.5;
% Construct H1
H1=H_MF_simple1(sig,omega,mu,eps,g);
[rho,timecell,fcell,dS]=find_bistab(H1,sig,J,gamma,Kappa,Rho_0,dT,rate,fmin,fmax,twait,resstep,extra_steps,method);
%%
title=['JCH bistability: J=' num2str(J) ' \gamma=\kappa=' num2str(gamma) '  \mu=' num2str(mu) '  \epsilon=\omega - J']
saveplot=true
savename='test2.png'
parforplot_a2f2color(rho,timecell,fcell,dS,sig,plottitle,saveplot,savename)
%% Some quick tries of previously unresolved scans
%% Check new find_bistab function
mu=0.4;
g=0.6;
fmax=0.5;
% Construct H1
H1=H_MF_simple1(sig,omega,mu,eps,g);
[rho,timecell,fcell,dS]=find_bistab(H1,sig,J,gamma,Kappa,Rho_0,dT,rate,fmin,fmax,twait,resstep,extra_steps,method);
%%
title=['JCH bistability: J=' num2str(J) ' \gamma=\kappa=' num2str(gamma) '  \mu=' num2str(mu) '  \epsilon=\omega - J']
saveplot=true
savename='test.png'
parforplot_a2f2color(rho,timecell,fcell,dS,sig,plottitle,saveplot,savename)
%% test SSP
mu=0.75;
g=0.5;
fmax=0.5;
savelocation='syms';
vplot=[];
[rho,timecell,fcell,dS]=SSP(mu,g,J,dT,rate,fmin,fmax,twait,resstep,extra_steps,method,vplot,savelocation,true,savelocation)
% test 2
mu=0.4;
g=0.6;
fmax=0.5;
savelocation='syms';
vplot=[];
[rho,timecell,fcell,dS]=SSP(mu,g,J,dT,rate,fmin,fmax,twait,resstep,extra_steps,method,vplot,savelocation,true,savelocation)
%% |<a*a>| plot 
mu=0.75;
g=0.5;
fmax=0.15;
savelocation='syms';
vplot=[0 1 0];
[rho,timecell,fcell,dS]=SSP(mu,g,J,dT,rate,fmin,fmax,twait,resstep,extra_steps,method,vplot,savelocation,true,savelocation,'',false)
%%
rho_full=[rho{1} rho{2} rho{3}];
time_full=[timecell{1} timecell{2} timecell{3}];
f_full=[fcell{1} fcell{2} fcell{3}];
n_tot=expectation_val(length(time_full),rho_full,sig{3});
%% Plots a
t=time_full;
figure,
subplot(3,1,2)
plot(t,n_tot,'color','black')
xlabel('time(s)');ylabel('|<a*a>|');%title('<a>');
subplot(3,1,3)
f_v=f_full;
plot(f_v.^2,n_tot,'color','black')
xlabel('f^2');ylabel('|<a*a>|');%title('<a>');
subplot(3,1,1)
plot(t,f_v,'color','black')
xlabel('t');ylabel('|f|^2');

sgtitle(['Trajectories bistability: J=' num2str(J) '    \gamma=\kappa=' num2str(gamma) ' \mu=0.75 g=0.5'])
%% Vergelijking met trajectories
% Het lijkt erop dat gamma en kappa misschien vergeten zijn in de
% trajectory calculations. We zullen hier eens kijken wat het resulaat is
% van gamma=0=kappa.
% DIT IS NIET WAAR  kappa en gamma zitten in H1 in plaats van H2
%% |<a*a>| plot 
mu=0.75;
g=0.5;
fmax=0.5;
savelocation='syms';
vplot=[];
[rho,timecell,fcell,dS]=SSP(mu,g,J,dT,rate,fmin,fmax,twait,resstep,extra_steps,method,vplot,savelocation,true,savelocation)