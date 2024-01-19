function [y]=sop_ris_sim(gambar_SR,gambar_RD,gambar_RE,gambar_JE,M,N,P,Q,K,nb1,nb2,Rs)
Nsim=4e6;
%% defining system parameters
unb1=2^(-nb1)*pi;
phi_err_1=2*unb1*(rand(M,Nsim,P)-.5);




unb2=2^(-nb2)*pi;
phi_err_2=2*unb2*(rand(N,Nsim,Q)-.5);





% %% path loss terms
% fc=2.4e9; %% 2.4 GHz
% c=3e8; %% speed of light
% dx=c/fc/2;
% dy=dx;
% Gt=10^(5/10);
% Gr=Gt;

%% path-loss definitions

%%% link via RIS


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

% landa_AE=gambar_RE*N;
% 
% landa_JE=gambar_JE*N;

%% rayleigh fading 
h_SL=1/sqrt(2).*(randn(M,Nsim,P)+1i*randn(M,Nsim,P));
h_LR=1/sqrt(2).*(randn(M,Nsim,P)+1i*randn(M,Nsim,P));
h_RL=1/sqrt(2).*(randn(N,Nsim,Q)+1i*randn(N,Nsim,Q));
h_LD=1/sqrt(2).*(randn(N,Nsim,Q)+1i*randn(N,Nsim,Q));
h_JL=1/sqrt(2).*(randn(N,Nsim,Q)+1i*randn(N,Nsim,Q));
h_LE=1/sqrt(2).*(randn(N,Nsim,K)+1i*randn(N,Nsim,K));

%% Rician fading 
% h_SL=sqrt(1/(2*(K_SL+1))).*(randn(M,K)+j*randn(M,K))+sqrt(K_SL/(K_SL+1));
% h_LR=sqrt(1/(2*(K_LR+1))).*(randn(M,K)+j*randn(M,K))+sqrt(K_LR/(K_LR+1));
% h_RL=sqrt(1/(2*(K_RL+1))).*(randn(N,K)+j*randn(N,K))+sqrt(K_RL/(K_RL+1));
% h_LD=sqrt(1/(2*(K_LD+1))).*(randn(N,K)+j*randn(N,K))+sqrt(K_LD/(K_LD+1));
% h_JL=sqrt(1/(2*(K_JL+1))).*(randn(N,K)+j*randn(N,K))+sqrt(K_JL/(K_JL+1));
% h_LE=sqrt(1/(2*(K_LE+1))).*(randn(N,K)+j*randn(N,K))+sqrt(K_LE/(K_LE+1));

gamma_SR=zeros(P,Nsim);
gamma_RD=zeros(Q,Nsim);
gamma_RE=zeros(K,Nsim);

for p=1:P
gamma_SR(p,:)=gambar_SR.*(abs(sum(abs(h_SL(:,:,p)).*abs(h_LR(:,:,p)).*exp(-j.*phi_err_1(:,:,p)),1))).^2;
end

for q=1:Q
gamma_RD(q,:)=gambar_RD.*(abs(sum(abs(h_RL(:,:,q)).*abs(h_LD(:,:,q)).*exp(-j.*phi_err_2(:,:,q)),1))).^2;
end

for k=1:K
gamma_RE(k,:)=gambar_RE.*(abs(sum(abs(h_RL(:,:,k)).*h_LE(:,:,k).*exp(-j.*(phi_err_2(:,:,k)+angle(h_LD(:,:,k)))),1))).^2./...
    (gambar_JE.*(abs(sum(h_JL(:,:,k).*h_LE(:,:,k).*exp(-j.*(phi_err_2(:,:,k)+angle(h_RL(:,:,k))+angle(h_LD(:,:,k)))),1))).^2+1);
end
gamma_SR_tot=sum(gamma_SR,1);
% sz_sr=size(gamma_SR_tot)
gamma_RD_tot=sum(gamma_RD,1);
% sz_rd=size(gamma_RD_tot)

if(K~=1)
gamma_RE_tot=max(gamma_RE);
else
gamma_RE_tot=gamma_RE;
end
% sz_re=size(gamma_RE_tot)

corrcoef(gamma_RD_tot,gamma_RE_tot)
y=sum(min(gamma_SR_tot,gamma_RD_tot)<(2^(2*Rs)*(1+gamma_RE_tot)-1))/Nsim;
