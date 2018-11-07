function phix  = phi_eval_1d(theta,xi,Llx,eta1,K,s,cg,Om,om,k0,sig,ep)
    Kvec = pi/Llx*[0:K -K+1:-1]';
    Xmesh = (-Llx:Llx/K:Llx-Llx/K)';
    Dx = 1i*Kvec;
    
    % Assume eta is in frequency space.
    eta1p = ifft(eta1);
    eta1sq = eta1p.^2;
    abseta1 = real(eta1p.*conj(eta1p));
    eta1xp = ifft(Dx.*eta1);
    theta2 = theta.^2;
    
    [a1,a2,a3,b3] = nls_expan_params(k0,Om,om,sig);
    
    n2 = a1*eta1sq;
    
    % non-quadrature based scheme
    
    zhoxv = -2*ep*Om*k0*abseta1;
    fhoxv = -s*Om*eta1p + 1i*ep*cg*eta1xp;
    shoxv = ep*(-2*Om*s*n2 + (Om -s*om)*k0*eta1sq);
    
    phixv = zhoxv + 2*real(fhoxv*theta + shoxv*theta2);        
    phix = spline(Xmesh,phixv,xi-Llx);
end

