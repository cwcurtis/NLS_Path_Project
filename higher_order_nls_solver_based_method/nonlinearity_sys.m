function rhs = nonlinearity_sys(w,k0,cg,om,Om,ep,sig,Dx,L0i,KT)

s = sign(k0);

eta0 = w(1:KT);
eta1 = w(KT+1:2*KT);

eta0p = real(ifft(eta0));
eta1p = ifft(eta1);

dxeta0p = real(ifft(Dx.*eta0));
Hdxeta0p = real(ifft(-abs(Dx).*eta0));

dxeta1p = ifft(Dx.*eta1);
meta1p = real(eta1p.*conj(eta1p));
Hdxmeta1p = real(ifft(-abs(Dx).*fft(meta1p)));

ellz = -k0*(om^2+2*Om^2-4*s*om*Om)/(2*(k0+om*Om-2*s*Om^2+4*sig*k0^3));
ello = (2*cg*k0*(Om-s*om)-2*s*om*Om+Om^2+om^2/2 + ellz*(1+om*cg+12*k0^2*sig-4*s*cg*Om))/(k0+om*Om-2*s*Om^2+4*sig*k0^3);

a0 = (k0^2*(-3/2*k0^3*sig+2*s*Om^2-s*om^2) + ellz*k0*(om^2+2*Om^2-4*s*om*Om))/(2*s*Om-om);
a1 = (k0*(4*cg*s*k0*Om-9*sig*k0^3+6*s*Om^2-3*s*om^2) + ellz*(2*om^2+4*Om^2-8*s*om*Om+4*k0*cg*(Om-s*om)) + ello*k0*(4*s*om*Om-2*Om^2-om^2))/(2*s*Om-om);
a2 = (k0*(-3/2*k0^3*sig+2*s*Om^2-s*om^2) + ellz*(om^2+2*Om^2-4*s*om*Om))/(2*s*Om-om);
a3 = k0*(2*s*Om-om);

fvec = a1*meta1p.*dxeta1p + a2*eta1p.^2.*conj(dxeta1p) + 1i*a3*eta1p.*Hdxmeta1p ... 
       + om*(om-2*s*Om-s*k0*cg)/(2*s*Om-om)*eta1.*dxeta0p ...
       + om*(om-2*s*Om-2*s*k0*cg)/(2*s*Om-om)*eta0p.*dxeta1p - 1i*cg*k0*eta1p.*Hdxeta0p;

zavec = om*(om-2*s*Om)*meta1p + ep*(cg*(om-2*s*Om)*Hdxmeta1p-2*s*om*cg*imag(conj(eta1p).*dxeta1p));   
   
rhs = [-Dx/(ep*om).*fft(zavec);L0i.*fft(1i*a0*meta1p.*eta1p - 1i*k0*om*eta0p.*eta1p + ep*fvec)];

%rhs = fft(1i*anl*meta1p.*eta1p);

    