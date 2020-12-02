%constants
gamma=0.005;
Kappa=0.005;
J=0.126;
f=@(t) 0.005;
w=1.11;
epsilon=w-J;
g=1;
mu=0;
fill=0;
Nx=4;
Ny=4;
Nxc=1;
Nyc=1;
T=2000;
dT=50;
spin=[0;1];
N_max=5;


n_sites=Nx*Ny;          %number of sites
Nc=Nxc*Nyc;             %number of sites in cluster
s_WF=2*(N_max+1);       %size wave function

sig = GetAllOperatorsSparse(Nc,N_max);
%%
tic
[time, Ct]=CalculateTrajectory(gamma, Kappa, J,f,w,epsilon,g,mu,Nxc,Nyc,Nx,Ny,T,dT,fill,5,spin);
toc
%% calculate averages
v_ps=time;      %vector containing the average photons/site at each time step
v_ss=time;      %vector containing the average spin/site at eacht time step
v_as=time;
for t = 1:length(time)
    [ps,ss,as] = GetAverageExpvals(Ct,n_sites,s_WF,t,sig);
    v_ps(t)=abs(ps);
    v_ss(t)=ss;
    v_as(t)=as;
end

%% plot averages
fig=figure,
%plot(time,v_ps,'.','MarkerSize',20,'Color','black')
plot(time,v_ps,'Color','black')
ylabel('Average amount of photons per site')
xlabel('Time step')
titlestr=['Avsps_f' num2str(f(0)) '_Nx' num2str(Nx) '_Ny' num2str(Ny)];
titlestr=[titlestr '_fill' num2str(fill) 'spin' num2str(spin(1)) '_' num2str(spin(2)) '.png'];
%saveas(fig,titlestr)

fig2=figure,
% plot(time,v_ss,'.','MarkerSize',20,'Color','black')
plot(time,v_ss,'Color','black')
ylabel('Average spin per site')
xlabel('Time step')
titlestr=['Avsss_f' num2str(f(0)) '_Nx' num2str(Nx) '_Ny' num2str(Ny)];
titlestr=[titlestr '_fill' num2str(fill) 'spin' num2str(spin(1)) '_' num2str(spin(2)) '.png'];
%saveas(fig2,titlestr)

fig3=figure,
% plot(time,v_as,'.','MarkerSize',20,'Color','black')
plot(time,v_as,'Color','black')
ylabel('Average |<a>| per site')
xlabel('Time step')
titlestr=['Avsas_f' num2str(f(0)) '_Nx' num2str(Nx) '_Ny' num2str(Ny)]
titlestr=[titlestr '_fill' num2str(fill) 'spin' num2str(spin(1)) '_' num2str(spin(2)) '.png']
%saveas(fig3,titlestr)





%%
% t=11;
% [ps,ss] = GetAverageExpvals(Ct,n_sites,s_WF,t,sig)
% ps=0;                   %average amount of photons per site
% ss=0;                   %average spin per site
% for i = 1:n_sites
%     Psi = Ct(t,(1+(i-1)*s_WF):(s_WF+(i-1)*s_WF))';      %pick wave function of site i. 
%                                                         %We transpose Psi for a more natural feeling(Psi is a ket)
%     ps = ps + Psi'*sig{1,3}*Psi;                        %expectation value counting operator
%     ss = ss + Psi'*sig{1,4}*Psi;           
% end
% ps=ps/n_sites
% ss=ss/n_sites