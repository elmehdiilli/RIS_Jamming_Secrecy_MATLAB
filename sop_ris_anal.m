function [yex,yan,yas]=sop_ris_anal(gambar_SR,gambar_RD,gambar_RE,gambar_JE,M,N,P,Q,K,nb1,nb2,Rs,n)
%                                 PNR_AB,JNR,v_AB,v_AE,dab,dar,drb,dae,dre,djr,djb,dje,dsj,dx,dy,Gt,Gr,bat_cap,Rs
%% defining fading parameters
unb1=2^(-nb1)*pi;
phi1_1=sin(unb1)/unb1;
phi1_2=sin(2*unb1)/(2*unb1);
% phi1_1=1;
% phi1_2=1;
mom_fad_SR=sqrt(pi/4);
mom_fad_RD=sqrt(pi/4);


m_SR=M/2*(phi1_1^2*mom_fad_SR^4)/(1+phi1_2-2*phi1_1^2*mom_fad_SR^4);
omega_SR=M^2*mom_fad_SR^4*phi1_1^2;


unb2=2^(-nb2)*pi;
phi2_1=sin(unb2)/unb2;
phi2_2=sin(2*unb2)/(2*unb2);
% phi2_1=1;
% phi2_2=1;
m_RD=N/2*(phi2_1^2*mom_fad_RD^4)/(1+phi2_2-2*phi2_1^2*mom_fad_RD^4);
omega_RD=N^2*mom_fad_RD^4*phi2_1^2;

% m_RD;


% %% path loss terms
% fc=2.4e9; %% 2.4 GHz
% c=3e8; %% speed of light
% dx=c/fc/2;
% dy=dx;
% Gt=10^(5/10);
% Gr=Gt;
% 
% %% path-loss definitions
% 
% %%% link via RIS
% 
% 
% teta_incid_S=abs(acos(x_S/d_SL));
% 
% teta_incid_R=abs(acos(x_R/d_RL));
% 
% 
% teta_incid_J=abs(acos(x_J/d_JL));
% 
% teta_refl_R=abs(acos(x_R/d_LR));
% teta_refl_D=abs(acos(x_D/d_LD));
% teta_refl_E=abs(acos(x_E/d_LE));
% 
% 
% PL_cst=Gt*Gr*dx^2*dy^2/(16*pi^2);
% 
% PL_RIS_SR=cos(teta_incid_S).*cos(teta_refl_R)*PL_cst./(d_SL*d_LR).^2;
% PL_RIS_RD=cos(teta_incid_R).*cos(teta_refl_D)*PL_cst./(d_RL*d_LD).^2;
% PL_RIS_RE=cos(teta_incid_R).*cos(teta_refl_E)*PL_cst./(d_RL*d_LE).^2;
% PL_RIS_JE=cos(teta_incid_J).*cos(teta_refl_E)*PL_cst./(d_JL*d_LE).^2;
% 
% %%% average SNR values
% gambar_SR=gambar_S*PL_RIS_SR;
% 
% gambar_RD=gambar_R*PL_RIS_RD;
% 
% gambar_RE=gambar_R*PL_RIS_RE;
% 
% gambar_JE=gambar_J*PL_RIS_JE;
% 
landa_RE=gambar_RE*N;
% 
landa_JE=gambar_JE*N;

%%%%%%%  A- evaluating the exact one (step 1)

ccdff=@(z) igamma(m_SR*P,m_SR/(omega_SR*gambar_SR)*z).*igamma(m_RD*Q,m_RD/(omega_RD*gambar_RD)*z)./(gamma(m_SR*P)*gamma(m_RD*Q));
yex=0;

for k=0:K-1
integrand=@(x) ccdff(2^(2*Rs)*(1+x)-1).*exp(-x/landa_RE).*(x./landa_RE+1/landa_JE+1)./...
(x./landa_RE+1/landa_JE).^2.*exp(-k.*x/landa_RE)./(x./landa_RE+1./landa_JE).^k;
yex=yex+nchoosek(K-1,k).*(-1)^k./(landa_JE^k).*integral(@(x) integrand(x),0,1e5);

end
yex=1-K/(landa_RE.*landa_JE).*yex;


% test PDF of illegitimate link:
% for k=0:K-1
% integrand=@(x) exp(-x/landa_RE).*(x./landa_RE+1/landa_JE+1)./...
% (x./landa_RE+1/landa_JE).^2.*exp(-k.*x/landa_RE)./(x./landa_RE+1./landa_JE).^k;
% yex=yex+nchoosek(K-1,k).*(-1)^k./(landa_JE^k).*integral(@(x) integrand(x),0,1e4);
% 
% end
% yex=K/(landa_RE.*landa_JE).*yex


%%%%%%%%%%%% B- Gauss-Laguerre form

%%% Gauss Laguerre parameters
% n=40;

syms t;
T=laguerreL(n,t);

Tp=sym2poly(T);

p=flip(roots(Tp));

ww=p./((n+1)^2*laguerreL(n+1,p).^2);


%%% Gauss Hermite parameters
% syms t;
% T=hermiteH(n,t);
% Tp=sym2poly(T);
% p=flip(roots(Tp));
% ww=2^(n-1)*sqrt(pi)*factorial(n)./(n^2*hermiteH(n-1,p).^2);

% ccdf_eq=@(z) (1-gammainc(m_SR*z/(omega_SR*gambar_SR),m_SR*P))*(1-gammainc(m_RD*z/(omega_RD*gambar_RD),m_RD*Q));
ccdf_eq=@(z) igamma(m_SR*P,m_SR*z/(omega_SR*gambar_SR))*igamma(m_RD*Q,m_RD*z/(omega_RD*gambar_RD));


sop=0;
sop_as=0;
for k=0:K-1
    for ii=1:n



% 
%   %% %% Gauss Laguerre : form 1
% 
%     func_integ=@(t) ccdf_eq(2^(2*Rs)*(landa_RE.*t+1)-1).*...
%     (t+1/landa_JE(1)+1)./ (t+1/landa_JE(1)).^(k+2).*exp(-k.*t);
%   sop=sop+ww(ii)*K/(landa_JE(1)^(k+1))*...
%     nchoosek(K-1,k)*(-1)^k*func_integ(abs(p(ii))); 
% 
% 
%       func_integ_as=@(t) ccdf_eq(2^(2*Rs)*(landa_RE.*t+1)-1).*...
%     (t+1/landa_JE(1)+1)./ (t+1/landa_JE(1)).^(k+2).*exp(-k.*t);
%   sop_as=sop_as+ww(ii)*K/(landa_JE(1)^(k+1))*...
%     nchoosek(K-1,k)*(-1)^k*func_integ(abs(p(ii))); 

% clc
  %% %% Gauss Laguerre : form 2
% 
%   func_integ=@(t) ccdf_eq(2^(2*Rs)*(landa_RE.*t/(k+1)+1)-1).*...
%     (t/(k+1)+1/landa_JE(1)+1)./ (t/(k+1)+1/landa_JE(1)).^(k+2);
%   sop=sop+ww(ii)*K/(landa_JE(1)^(k+1)*(k+1))*...
%     nchoosek(K-1,k)*(-1)^k*func_integ(p(ii)); 


  %% %% Gauss Laguerre : form 3
% 
    func_integ=@(t) ccdf_eq(2^(2*Rs)*(t+1)-1).*...
    (t/landa_RE+1/landa_JE(1)+1)./ (t/landa_RE+1/landa_JE(1)).^(k+2).*exp(-(k+1).*t/landa_RE).*exp(t);
  sop=sop+ww(ii)*K/(landa_RE*landa_JE(1)^(k+1))*...
    nchoosek(K-1,k)*(-1)^k*func_integ(p(ii)); 

      func_integ_as=@(t) (m_SR/omega_SR*(2^(2*Rs)*(t+1)-1)).^(m_SR*P)/gamma(m_SR*P+1).*...
    (t/landa_RE+1/landa_JE(1)+1)./ (t/landa_RE+1/landa_JE(1)).^(k+2).*exp(-(k+1).*t/landa_RE).*exp(t);
  sop_as=sop_as+ww(ii)*K/(landa_RE*landa_JE(1)^(k+1))*...
    nchoosek(K-1,k)*(-1)^k*func_integ_as(p(ii)); 
% 
% 

%  %% Gauss Hermite
% % 
%   func_integ=@(z) ccdf_eq(2^(2*Rs)*(landa_RE.*exp(z)+1)-1).*...
%     (exp(z)/landa_RE+1/landa_JE(1)+1)./ (exp(z)/landa_RE+1/landa_JE(1)).^(k+2).*exp(z).*exp(-(k+1).*exp(z)/landa_RE).*exp(z.^2);
%   sop=sop+ww(ii)*K/(landa_JE(1)^(k+1)*landa_RE)*...
%     nchoosek(K-1,k)*(-1)^k*func_integ(p(ii)); 



    end
end

yan=1-sop/(gamma(m_SR*P)*gamma(m_RD*Q));
yas=sop_as*gambar_SR^(-m_SR*P);

% yan=1-sop;

%%%%%%%  C- evaluating the triple complex integr. (N.B: it was implemented
%%%%%%%  for the previous scheme -- SISO !!!!!)

% term11=0;
% term12=0;
% 

% term21=0;
% term22=0;
% for k=0:lmax
%     for l=0:lmax
% 
% f11=@(s,v) gammaz(s).*gammaz(m_SR+s).*gammaz(v).*gammaz(m_RD+v)./(gammaz(1+s).*gammaz(1+v))...
%     .*gammaz(-s-v-l).*gammaz(s+v+k).*gammaz(s+v+l+k+1)./(gammaz(s+v).*gammaz(s+v+k+1)).*...
%     (N*gambar_JE*m_SR*(2^(2*Rs)-1)/(omega_SR*gambar_SR)).^(-s).*(N*gambar_JE*m_RD*(2^(2*Rs)-1)/(omega_RD*gambar_RD)).^(-v);
% 
% f12=@(s,v) gammaz(s).*gammaz(m_SR+s).*gammaz(v).*gammaz(m_RD+v)./(gammaz(1+s).*gammaz(1+v))...
%     .*gammaz(s+v-l).*gammaz(s+v+k).*gammaz(k+l+1)./(gammaz(s+v).*gammaz(s+v+k+1)).*...
%     (m_SR*(2^(2*Rs)-1)/(omega_SR*gambar_SR)).^(-s).*(m_RD*(2^(2*Rs)-1)/(omega_RD*gambar_RD)).^(-v);
% 
% 
% f21=@(s,v) gammaz(s).*gammaz(m_SR+s).*gammaz(v).*gammaz(m_RD+v)./(gammaz(1+s).*gammaz(1+v))...
%     .*gammaz(-s-v-l-1).*gammaz(s+v+k).*gammaz(s+v+l+k+2)./(gammaz(s+v).*gammaz(s+v+k+2)).*...
%     (N*gambar_JE*m_SR*(2^(2*Rs)-1)/(omega_SR*gambar_SR)).^(-s).*(N*gambar_JE*m_RD*(2^(2*Rs)-1)/(omega_RD*gambar_RD)).^(-v);
% 
% f22=@(s,v) gammaz(s).*gammaz(m_SR+s).*gammaz(v).*gammaz(m_RD+v)./(gammaz(1+s).*gammaz(1+v))...
%     .*gammaz(s+v+1-l).*gammaz(s+v+k).*gammaz(k+l+1)./(gammaz(s+v).*gammaz(s+v+k+2)).*...
%     (m_SR*(2^(2*Rs)-1)/(omega_SR*gambar_SR)).^(-s).*(m_RD*(2^(2*Rs)-1)/(omega_RD*gambar_RD)).^(-v);
% 
% 
% cs=1/2;
% cv=1/2;
% 
% term11=term11+1/gambar_JE*(-1)^l/(factorial(l)*factorial(k))*(1/(N*gambar_JE))^l*...
%     integral(f11,cs-20*i,cs+20*i,cv-20*i,cv+20*i)*(1-gambar_RE/(gambar_JE*(2^(2*Rs)-1)/2^(2*Rs)))^k;
% 
% term12=term12+1/gambar_JE*(-1)^l/(factorial(l)*factorial(k))*(1/(N*gambar_JE))^l*...
%     integral(f12,cs-20*i,cs+20*i,cv-20*i,cv+20*i)*(1-gambar_RE/(gambar_JE*(2^(2*Rs)-1)/2^(2*Rs)))^k;
% 
% term21=term21+N*(-1)^l/(factorial(l)*factorial(k))*(1/(N*gambar_JE))^(l+1)*...
%     integral(f21,cs-20*i,cs+20*i,cv-20*i,cv+20*i)*(1-gambar_RE/(gambar_JE)*(2^(2*Rs)-1)/2^(2*Rs))^k;
% 
% term22=term22+N*(-1)^l/(factorial(l)*factorial(k))*(1/(N*gambar_JE))^l*...
%     integral(f22,cs-20*i,cs+20*i,cv-20*i,cv+20*i)*(1-gambar_RE/(gambar_JE)*(2^(2*Rs)-1)/2^(2*Rs))^k;
%     end
% end
% 
% yan=1-1/((2*pi*j)^3*N*gamma(m_SR)*gamma(m_RD))*(term11+term12+term21+term22);


end