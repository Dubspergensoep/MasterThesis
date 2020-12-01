function H=H_MF_simple1(sig,w,mu,eps,g)
    %time independent
    H = (w - mu)*sig{3} + 0.5*(eps-mu)*sig{4} ...
          + g*(sig{2}*sig{5} + sig{1}*sig{6});
end