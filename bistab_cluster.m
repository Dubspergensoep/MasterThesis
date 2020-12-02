%% Bistability cluster
%% Constants
% system/hamiltonian parameters
gamma=0.005;
Kappa=0.005;
J=0.126;
omega=1;
epsilon=omega-J;
mu=0.75;
g=0.5;

% sym parameters
fill=0;
Nx=4;
Ny=4;
Nxc=1;
Nyc=1;
spin=[0;1];
N_max=5;
n_sites=Nx*Ny;          %number of sites
Nc=Nxc*Nyc;             %number of sites in cluster
s_WF=2*(N_max+1);       %size wave function
sig = GetAllOperatorsSparse(Nc,N_max);
% Laser parameters
fmax=0.15;
rate1=1e-5;
rate2=1e-5;
twait=1000;
t_fmax=fmax/rate1;
f=@(t) laser(t,rate1,t_fmax,twait,rate2)
t_end=t_fmax+twait+fmax/rate2;
dT=10;

%% Calculate cluster
parpool(4)
tic
parfor j=1:28 %Since there are 28 notes we will use all of them
    [time, Ct]=CalculateTrajectory(gamma, Kappa, J,f,omega,epsilon,g,mu,Nxc,Nyc,Nx,Ny,t_end,dT,fill,5,spin);
    %Calculate |<a>|
    v_ea=[];
    for n=1:length(time)
        ea=0;                                               %expectation value a
        for i = 1:n_sites
            Psi = Ct(n,(1+(i-1)*s_WF):(s_WF+(i-1)*s_WF))';      %pick wave function of site i.
            ea = ea + abs(Psi'*sig{1,2}*Psi);
        end
        v_ea(n)=abs(ea)/n_sites;
    end
    fname=strjoin(["sysm/syms_" num2str(j)],'');
    parsave(fname,v_ea,time,Ct)
end
toc

%% functions
function parsave(fname, v_ea,time,Ct)
  save(fname, 'v_ea', 'time','Ct')
end