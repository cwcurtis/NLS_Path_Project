function [xtrack,ztrack,w] = nls_solver(K,ad,anl,cg,k0,kval,Om,om,sig,ep,tf,dt,flag)

    % solves u_t = dk2Om*u_xx + a3*|u|^2*u 

    % Llx is domain size, i.e. solve on domain [-Llx,Llx]
    % sig is -1 for bright/focusing, +1 for dark/defocusing
    % tf is total simulation time
    % mu is magnitude of damping on envelope of dark initial condition
        KT = 2*K;
        nmax = round(tf/dt); % Step size for time integrator and number of time steps
    
        Tval = 2*ellipke(kval^2);
        X = (-Tval:Tval/K:Tval-Tval/K)'; % Spatial mesh
        Kvec = pi/Tval*[0:K -K+1:-1]';
        Dx = 1i*Kvec;
        Dx2 = Dx.^2;
        scfac2 = dt*pi/Tval;
        cgsc = pi*cg/(ep*Tval);
        s = sign(k0);
        L0i = (1 - 2*1i*ep*s*cg/(om-2*s*Om)*Dx);
        Lap = dt*(L0i*1i*ad+ ep*sig/(om-2*s*Om)*Dx).*Dx2;
        avspd = om^2*(2*s*Om-om)/(1+cg*om) - 2*k0*Om*(1 + (1-s*om/Om)^2);
        % Bright/Focusing initial conditions

        if flag == 1
            if ad/anl > 0
                [~,cnvals,~] = ellipj(X,kval^2);
                uint = kval*sqrt(2*ad/anl)*cnvals;             
            else
                [snvals,~,~] = ellipj(X,kval^2);
                uint = kval*sqrt(-2*ad/anl)*snvals;             
            end
            disp('Lagrangian Speed')
            disp(avspd*2*abs(ad/anl))
            disp('Stokes Drift Speed')
            disp((-2*k0*Om*(1 + (1-s*om/Om)^2))*(2*abs(ad/anl)))        
        else
            uint = ones(KT,1);
            disp('Lagrangian Speed')
            disp(avspd*2)
            disp('Stokes Drift Speed')
            disp((-2*k0*Om*(1 + (1-s*om/Om)^2))*2)        
        end
        
        [a1,a2,a3,b3] = nls_expan_params(k0,Om,om,sig);
        n00 = om*(2*Om*s-om)/(1+cg*om);
                
        uf = fft(uint); 
        
        n0 = n00*abs(uint).^2;
        n2 = a1*uint.^2;
        sint = ep^2*n0 + 2*ep*real(uint.*exp(1i*k0*X/ep) + ep*n2.*exp(2*1i*k0*X/ep));           
            
        Nvorts = 2;
        Lind = K;
        Rind = K+2;
        
        xint1 = pi/Tval*(X(Lind) + Tval);
        xint2 = pi/Tval*(X(Rind) + Tval);
        
        xpos = [xint1; xint2];
        zpos = [sint(Lind);sint(Rind)];
        
        xtrack = zeros(Nvorts,nmax+1);
        ztrack = zeros(Nvorts,nmax+1);
        
        xtrack(:,1) = (-Tval + Tval/pi*xpos)/ep;
        ztrack(:,1) = zpos;       
        
        w = uf;
        
    % Begin setup of RK4 method

        E = exp(Lap/2);
        E2 = E.^2;
                
    % Solve NLS equation in time
    
        for nn=1:nmax
            
            t = (nn-1)*dt;
            
            wp = ifft(w);
            a = dt*nonlinearity(w,k0,cg,om,Om,ep,sig,ad,anl,Dx);
            af = E.*(w+a/2);
            
            wap = ifft(af);
            b = dt*nonlinearity(af,k0,cg,om,Om,ep,sig,ad,anl,Dx);
            bf = E.*w + b/2;
            
            wbp = ifft(bf);
            c = dt*nonlinearity(bf,k0,cg,om,Om,ep,sig,ad,anl,Dx);
            cf = E2.*w + E.*c;
            
            wcp = ifft(cf);
            d = dt*nonlinearity(cf,k0,cg,om,Om,ep,sig,ad,anl,Dx);
            
            for jj = 1:Nvorts
            
                xwave = Tval/pi*xpos(jj) - Tval;
                xlocpos = xloc_calc(xpos(jj),cgsc,t,Tval);
                theta = exp(1i*(k0*xwave + Om*t/ep)/ep);
                indc = loc_find(xlocpos,X);
                wval = wp(indc);
                
                n0 = n00*abs(wval.^2);
                n2 = a1*wval.^2;
                zval = ep^2*n0 + 2*ep*real(wval.*exp(1i*(k0*xwave/ep+ Om*t/ep^2)) + ep*n2.*exp(2*1i*(k0*xwave/ep+ Om*t/ep^2)) );
                
                xdot = phi_eval_1d(theta,Tval*(xpos(jj)+cgsc*t)/pi,Tval,w,K,s,cg,Om,om,k0,sig,ep);
                
                ax = scfac2*(om*zval/ep + xdot);
                
                xwave = Tval/pi*(xpos(jj)+ax/2) - Tval;
                xlocpos = xloc_calc(xpos(jj)+ax/2,cgsc,t+dt/2,Tval);
                theta = exp(1i*(k0*xwave + Om*(t+dt/2)/ep)/ep);
                indc = loc_find(xlocpos,X);
                wval = wap(indc);
                
                n0 = n00*abs(wval.^2);
                n2 = a1*wval.^2;
                zval = ep^2*n0 + 2*ep*real(wval.*exp(1i*(k0*xwave/ep+ Om*t/ep^2)) + ep*n2.*exp(2*1i*(k0*xwave/ep+ Om*t/ep^2)) );
                            
                xdot = phi_eval_1d(theta,Tval*(xpos(jj)+cgsc*(t+dt/2)+ax/2)/pi,Tval,af,K,s,cg,Om,om,k0,sig,ep);
                 
                bx = scfac2*(om*zval/ep + xdot);
                
                xwave = Tval/pi*(xpos(jj) + bx/2) - Tval;
                xlocpos = xloc_calc(xpos(jj)+bx/2,cgsc,t+dt/2,Tval);
                theta = exp(1i*(k0*xwave + Om*(t+dt/2)/ep)/ep);
                indc = loc_find(xlocpos,X);
                wval = wbp(indc);
                
                n0 = n00*abs(wval.^2);
                n2 = a1*wval.^2;
                zval = ep^2*n0 + 2*ep*real(wval.*exp(1i*(k0*xwave/ep+ Om*t/ep^2)) + ep*n2.*exp(2*1i*(k0*xwave/ep+ Om*t/ep^2)) );
                
                xdot = phi_eval_1d(theta,Tval*(xpos(jj)+cgsc*(t+dt/2)+bx/2)/pi,Tval,bf,K,s,cg,Om,om,k0,sig,ep);
                                
                cx = scfac2*(om*zval/ep + xdot);
                
                xwave = Tval/pi*(xpos(jj) + cx) - Tval;
                xlocpos = xloc_calc(xpos(jj)+cx,cgsc,t+dt,Tval);
                theta = exp(1i*(k0*xwave + Om*(t+dt)/ep)/ep);
                indc = loc_find(xlocpos,X);
                wval = wcp(indc);
                
                n0 = n00*abs(wval.^2);
                n2 = a1*wval.^2;
                zval = ep^2*n0 + 2*ep*real(wval.*exp(1i*(k0*xwave/ep+ Om*t/ep^2)) + ep*n2.*exp(2*1i*(k0*xwave/ep+ Om*t/ep^2)) );
                
                xdot = phi_eval_1d(theta,Tval*(xpos(jj)+cgsc*(t+dt)+cx)/pi,Tval,cf,K,s,cg,Om,om,k0,sig,ep);
                
                dx = scfac2*(om*zval/ep + xdot);
                
                xpos(jj) = xpos(jj) + (ax + 2*(bx+cx) + dx)/6;
                zpos(jj) = zval;                  
                
                
            end    
            
            w = E2.*w + (E2.*a + 2*E.*(b+c) + d)/6;                
            
            xtrack(:,nn+1) = (-Tval + Tval*xpos/pi)/ep;
            ztrack(:,nn+1) = zpos;            
            
        end
        
        wp = ifft(w);
        no_cells = 1+mod(round(abs(cg)*(t+dt)/ep*K/Tval),2*K);
        wp = wp([no_cells:2*K 1:no_cells-1]');
        n0p = n0*abs(wp.^2);
        n2 = a1*wp.^2;
        w = ep^2*n0p + 2*ep*real( wp.*exp(1i*(k0*X/ep+ Om*(t+dt)/ep^2)) + ep*n2.*exp(2*1i*(k0*X/ep+ Om*(t+dt)/ep^2)) );           
        
end