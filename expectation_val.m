function v_ea=expectation_val(l_time,rho_l_data,A)
    s_WF=sqrt(length(rho_l_data(:,1)));
    v_ea=zeros(1,l_time);
    for i=1:l_time
        rho = reshape(rho_l_data(:,i), s_WF, s_WF);
        v_ea(i)=abs(trace(A*rho));
    end
end