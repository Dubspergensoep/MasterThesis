function [ps,ss]=GetAverageExpvals(Ct,n_sites,s_WF,t,sig)
%Ct result - from CalculateTrajectory functions containing all wave functions
%n_sites   - amount of sites
%s_WF      - size single site wave function
%t         - timestep at which the average should be calculated
%sig       - Containts operators technacally this could be calculated here

ps=0;                   %average amount of photons per site
ss=0;                   %average spin per site
for i = 1:n_sites
    Psi = Ct(t,(1+(i-1)*s_WF):(s_WF+(i-1)*s_WF))';      %pick wave function of site i. 
                                                        %We transpose Psi for a more natural feeling(Psi is a ket)
    ps = ps + Psi'*sig{1,3}*Psi;                        %expectation value counting operator
    ss = ss + Psi'*sig{1,4}*Psi;           
end
ps=ps/n_sites;
ss=ss/n_sites;
end