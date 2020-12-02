%recreate figures
%recreation figure 5.5
diary remove_1_gamma
gamma=0.005;
Kappa=0.005;
f=@(t) 0.005;
fill=0;
Nx=4;
Ny=4;
Nxc=1;
Nyc=1;
T=2000;                   %T=10 lijkt goed genoeg zelfde resultaat voor 20 &40
dT=50;
spin=[0;1];
N_max=5;

n_sites=Nx*Ny;          %number of sites
Nc=Nxc*Nyc;             %number of sites in cluster
s_WF=2*(N_max+1);       %size wave function
t=T+1;                  %we use this to select the last Psi from Ct so the initial condition does not matter asmuch
sig = GetAllOperatorsSparse(Nc,N_max);

g=1;
J=0.126*g;
%J=0.004*g;

mu=0;
%x=1;            %our variable (mu-omega)/g
% v_x=[-1.5:0.1:-1.1 -1.05:0.05:-0.9 -0.8:0.1:0.8 0.9:0.05:1.1 1.2:0.1:1.5];
% v_x=[-1.5:0.1:-0.9 -0.8:0.2:0.6 0.7:0.1:1.5];
% v_x=[-1.5:0.05:-1.25 -1.24:0.01:-1.05 -1:0.1:0.75 0.8:0.01:0.95 1:0.05:1.5];
v_x=[-1.5:0.1:-1.25 -1.24:0.01:-1.05 -1:0.4:0.75 0.8:0.01:0.95 1:0.1:1.5];
v_ea=v_x;

%% run all tests
parfor j=1:length(v_x)
    w=-(v_x(j)*g-mu);
    epsilon=w-J;
    [time, Ct]=CalculateTrajectory(gamma, Kappa, J,f,w,epsilon,g,mu,Nxc,Nyc,Nx,Ny,T,dT,fill,5,spin);
    %Calculate |<a>|
    v_eat=zeros(1,T/(2*dT)); %looking for a better name
    for k=1:length(v_eat)
        ea=0;
        for i = 1:n_sites
            Psi = Ct(T/(2*dT)+k,(1+(i-1)*s_WF):(s_WF+(i-1)*s_WF))';      %pick wave function of site i.
            ea = ea + abs(Psi'*sig{1,2}*Psi);    
        end
        v_eat(k)=abs(ea)/n_sites;
    end
    %disp(v_eat)
    v_ea(j)=mean(v_eat);
end
%% Make plots
figure
plot(v_x,v_ea,'.','MarkerSize',20,'Color','black')
xlabel('(\mu-\omega)/g')
ylabel('|<a>|')
title('figure 5.5 Nissen bottom left - with trajectories')
parameters=struct;
parameters.gamma=gamma;parameters.Kappa=Kappa;parameters.f=f;parameters.fill=fill;parameters.Nx=Nx;parameters.Ny=Ny;parameters.T=T;parameters.spin=spin;parameters.N_max=N_max;parameters.J=J;parameters.g=g;parameters.mu=mu
%%
save('sym_bottom_left_10x10_verplaats_remove1_gamma','parameters','v_ea','v_x')
diary off
%% full phase plot
v_J=0:0.01:0.15;
[m_J, m_x]= meshgrid(v_J, v_x);
m_ea=zeros(size(m_J));

for j=1:length(v_J)
    J=v_J(j);
    for c=1:length(v_x)
        w=-(v_x(c)*g-mu);
        epsilon=w-J;
        [time, Ct]=CalculateTrajectory(gamma, Kappa, J,f,w,epsilon,g,mu,Nxc,Nyc,Nx,Ny,T,1,fill,5,spin);
        %Calculate |<a>|
        ea=0;                                               %expectation value a
        for i = 1:n_sites
            Psi = Ct(t,(1+(i-1)*s_WF):(s_WF+(i-1)*s_WF))';      %pick wave function of site i.
            ea = ea + Psi'*sig{1,2}*Psi;        
        end
        m_ea(c,j)=abs(ea)/n_sites;  %store values
    end
end
%% Make phase plot
figure,
contourf(m_J,m_x,m_ea)
colorbar
xlabel('J/g')
ylabel('(\mu-\omega)/g')
map=[1 1 1
    0 0 0.9
    0 0 0.8
    0 0 0.7
    0 0 0.6
    0 0 0.5
    0 0 0.4]
colormap(map)



