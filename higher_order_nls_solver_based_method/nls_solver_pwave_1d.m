function [Tvals,xtrack,ztrack,w,sdrift] = nls_solver_pwave_1d(K,Llx,A,cg,anl,k0,Om,om,sig,ep,tf)
     
    X = (-Llx:Llx/K:Llx-Llx/K)'; % Spatial mesh
    s = sign(k0);
    n0 = om*(2*Om*s-om)/(1 + cg*om);
    disp(om*n0)
    [a1,a2,a3,b3] = nls_expan_params(k0,Om,om,sig);
    a31 = -( (3*k0*om^2+12*Om^2*k0-18*Om*k0*s*om)*a1 + 9/2*k0^5*sig - 3*Om^2*k0^2*s - 3*k0^2*s*om^2 + 6*Om*k0^2*om)/(3*k0 + 3*Om*om - 9*s*Om^2 + 27*sig*k0^3);
    
    Lind = K;
    Rind = K+2;
      
    xint1 = X(Lind);
    xint2 = X(Rind);
        
    xpos = [xint1; xint2];
    
    avspd = om^2*(2*s*Om-om)/(1+cg*om) - 2*k0*Om*(1 + (1-s*om/Om)^2);
    disp('Lagrangian Speed')
    disp(avspd*A^2)
    disp('Stokes Drift Speed')
    disp((-2*k0*Om*(1 + (1-s*om/Om)^2))*(A^2))
    
    % Solve NLS equation in time
    opts = odeset('RelTol',1e-11,'AbsTol',1e-12);
    [Tvals,xtrack] = ode45(@(t,lhs) phi_eval_pwave_ho_1d(t,lhs,s,A,anl,Om,om,n0,k0,sig,cg,ep),[0 tf],xpos(2),opts);
    %[Tvals,theta_track] = ode45(@(t,lhs) phi_eval_pwave_ho_1d_theta(t,lhs,s,A,anl,Om,om,n0,a1,k0,ep),[0 tf],k0*xpos(2),opts);
    %xtrack = (theta_track - (Om + anl*A^2*ep^2)*Tvals)/k0;
    
    %plot(Tvals,xtrack)
    %pause
    
    sdrift = zeros(length(Tvals),1);
    ztrack = zeros(length(Tvals),1);
    
    for jj=1:length(Tvals)
        tloc = Tvals(jj);
        xloc = xtrack(jj);
        w = A*exp(1i*anl*ep^2*A^2*tloc);    
        n2 = a1*w^2;
        n3 = a31*w^3;
        theta = exp(1i*(k0*xloc+Om*tloc));
        
        ztrack(jj) = ep^2*n0*A^2 + 2*ep*real( w*theta + ep*n2*theta^2 + ep^2*n3*theta^3);
        sdrift(jj) = Stokes_Drift_pwave(w,sign(k0),Om,om,k0);                             
    end
    
    w = A*exp(1i*anl*ep^2*A^2*tf);
    n2 = a1*w.^2;
    w = ep^2*n0*A^2 + 2*ep*real( w.*exp(1i*(k0*X + Om*tf)) + ep*n2.*exp(2*1i*(k0*X+ Om*tf)) );             
    
end