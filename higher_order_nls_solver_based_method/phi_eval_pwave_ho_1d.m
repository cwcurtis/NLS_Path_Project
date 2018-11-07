function rhs = phi_eval_pwave_ho_1d(t,lhs,s,A,anl,Om,om,n0,k0,sig,cg,ep)
    
    x = lhs;
    
    % Assume eta is in physical space.
    theta = exp(1i*(k0*x + Om*t));
    eta1p = A*exp(1i*A^2*anl*ep^2*t);
    
    eta1sq = eta1p.^2;
    eta1cu = eta1sq.*eta1p;          
    theta2 = theta.^2;
    theta3 = theta2.*theta;
    
    [a1,a2,a3,b3] = nls_expan_params(k0,Om,om,sig);
    a31 = -( (3*k0*om^2+12*Om^2*k0-18*Om*k0*s*om)*a1 + 9/2*k0^5*sig - 3*Om^2*k0^2*s - 3*k0^2*s*om^2 + 6*Om*k0^2*om)/(3*k0 + 3*Om*om - 9*s*Om^2 + 27*sig*k0^3);
    
    n2 = a1*eta1sq;
    n3 = a31*eta1cu;
    
    z = ep^2*n0*A^2 + 2*ep*real( eta1p*theta + ep*n2*theta2 + ep^2*n3*theta3);
    ztrun = 2*ep*real( eta1p*theta );      
    fhox = -s*Om*eta1p-ep^2/2*(k0^2*(Om*s-2*om)+ 2*k0*om*(s*om-Om)*(2*s*Om-om)/(1+om*cg) + 2*k0*a1*(s*om-3*Om) + 2*s*anl)*eta1p.*A^2;
    shox = ep*(-2*s*Om*a1+(2*Om-s*om)*k0)*eta1sq;
    thox = -ep^2/2*(6*Om*s*a31 + (6*k0*s*om-18*Om*k0)*a1 + (9*Om*k0^2*s-6*k0^2*om))*eta1cu;
     % Compute phix
    
    At11 = fhox.*theta;
    At22 = shox*theta2;  
    At33 = thox*theta3;
    
    Rk1 = (1 + abs(k0)*ztrun).*At11;
    Rk2 = At22;
    Rk3 = At33;
    phix = 2*real(Rk1 + Rk2 + Rk3);
    
    rhs = om*z + ep*phix;
    
end

