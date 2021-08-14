function [dBdt,struct1,struct2]= D_G_ode_func(t,B,D,G,T,seasonal) 

%This ODE function calculate all the variables (uptake rates, division
%rates, mortality) needed for the resolution of the differential equations 


global Imax Jmax x xo y 
global alphaF palat mp0 m2 Z 
global gmax P0 reminN reminS D0
global JR betaL betaN betaSi betaF
global L N0 S0 dil 

%==== Reshaping biomass into a matrix for calculation ====%
 disp(t);
 Bd=reshape(B(1:Jmax*Imax),[Jmax Imax]); %Diatoms
 D.b = Bd;
 Bg = reshape(B(Jmax*Imax+1:Jmax*2*Imax),[Jmax Imax]); %Generalists
 G.b = Bg;

 %==== Environmental conditions ====%
 if seasonal==1     %Light, dilution rate and grazers abundance for each day 
     N = B(2*Jmax*Imax+1);
     Si = B(2*Jmax*Imax+2);
     if roundn(t(end),0)==0
         Lday=L(1);
         dilday=dil(1);
         Zday = Z(1);
     else
         Lday=L(roundn(t(end),0));
         dilday=dil(roundn(t(end),0));
         Zday = Z(roundn(t(end),0));
     end
 else %for steady state
     Lday = L;
     dilday = dil;
     Zday = Z;
     N=N0;
     Si = S0;
 end 
 
%==== Trophic interactions ====%
 Fd = zeros(1,Imax); Fg = zeros(1,Imax);
 G.JF= zeros(Jmax,Imax);
 D.m = zeros(Jmax,Imax);     
 Bpred = sum(G.b.*G.phif);   %the biomass of predators [mg C / m-3] tin size classes r(i)
 for i = 1:Imax  % calculated once
    Fd(i) = sum(T(:,:,i).*D.b,'all');%food available from diatoms [mg C / m-3] to all generalists in size class x(i)
    Fg(i) = sum(T(1,:,i).*G.b,'all'); %food available from generalists [mg C / m-3] to all generalists in size class x(i)
    G.JF(:,i) = alphaF*G.phif(:,i)*(palat*Fd(i) + Fg(i));% J_F=A_F*F in size class x(i)
    %D.m = D.m + alphaF*Bpred(i)*T(:,:,i); %predation by larger cells for diatoms 
 end

%   G.m = D.m(1,:);
%   D.m = palat*D.m; % mortality rates [day-1] due to internal predation 
%   D.m(D.m<0)=0;G.m(G.m<0)=0;

  %==== Uptake rates ====%
  G.JL = Lday.*G.AL; G.JN=N*G.AN;
    
  D.JL = Lday.*D.AL;  D.JN=N*D.AN; D.JSi= Si*D.AN;

 %==== Downregulation ====%
 D.Jresp = JR.*x; %Basal metabolism
 G.Jresp = JR.*x;
 
 D.limind = ones(Jmax,Imax);
 D.epsL= ones(Jmax,Imax);D.epsN = zeros(Jmax,Imax);D.epsSi = D.epsN;
 
 G.eps = max(0,min(1,(G.JL*(1-betaL)-G.Jresp-betaF*G.JF)./((betaN+G.CN).*G.JN)));

 %Different eventual cases :
 D.epsN1 = max(0,min(((1-betaL)*D.JL-D.Jresp)./((betaN + D.CN.*(1+betaSi./D.CS)).*D.JN),1)); %C limiting
 D.epsSi1 = max(0,min(((1-betaL)*D.JL-D.Jresp)./((betaSi + D.CS + betaN.*D.CS./D.CN).*D.JSi),1));
    
 D.epsSi2 = max(0,min(D.CN.*D.JN./(D.JSi.*D.CS),1)); %N limiting 
 D.epsL2 = max(0,(D.Jresp + D.CN.*D.JN +betaN.*D.JN + betaSi*D.CN.*D.JN./D.CS)./((1-betaL)*D.JL));
    
 D.epsN3 = max(0,min(1,D.CS.*D.JSi./(D.CN.*D.JN)));%Si limiting
 D.epsL3 =max(0,(D.Jresp + (D.CS+betaSi).*D.JSi + betaN*D.CS.*D.JSi./D.CN)./((1-betaL)*D.JL));    
    
 for i=1:Jmax*Imax
    if D.epsN1(i)<1 && D.epsSi1(i)<1 %carbon is limiting
        D.epsN(i) = D.epsN1(i);
        D.epsSi(i) = D.epsSi1(i);
    elseif D.epsSi2(i)<1 && D.epsL2(i) <1 %Nitrogen is limiting
        D.epsSi(i) = D.epsSi2(i); 
        D.limind(i) = 2 ; 
        D.epsN(i) = 1; 
    elseif D.epsN3(i)<1 && D.epsL3(i)<1 %Silica is limiting 
        D.epsN(i) = D.epsN3(i);
        D.epsSi(i) = 1; 
        D.limind(i) = 3;      
    end           
 end 
     
%==== Effective uptakes ====%
D.JNreg = D.epsN.*D.JN;
D.JSireg = D.epsSi.*D.JSi;
D.Ccons = max(0,(1-betaL)*D.JL-D.Jresp-betaL*D.epsL.*D.JL-betaN*D.JNreg-betaSi*D.JSireg); 
D.Ncons = D.CN.*D.JNreg;
D.Sicons = D.CS.*D.JSireg;
D.Jeff = min(min(D.Ncons,D.Sicons),D.Ccons); 
D.Jeff = max(D.Jeff,0); %Effective uptake (cannot be negative)

    
G.JNreg = G.JN.*G.eps ; 
G.Ccons = max(0,G.JF+(1-betaL)*G.JL-G.Jresp-betaF*G.JF-betaN*G.JNreg);
G.Ncons = G.CN.*G.JNreg+G.JF ;
G.Jeff = min(G.Ccons,G.Ncons);

%==== Cells division rate ====%
D.g = gmax.*D.Jeff./(D.Jeff+gmax.*ones(Jmax,Imax));%division rate
G.g = gmax.*G.Jeff./(G.Jeff+gmax.*ones(Jmax,Imax));

%==== Grazing by copepods ====%
 D.mp = zeros(Jmax,Imax); G.mp = zeros(Jmax,Imax);
    for i=1:Jmax*Imax
        if x(i)> 10^(-2)
            D.mp(i)= mp0*(1-palat)*Zday*D.r(i).^(-3/4);
            G.mp(i)= mp0*Zday*G.r(i).^(-3/4);
        end 
    end 

%==== Effective uptake fluxes of N and Si ====%
G.etaN = zeros(Jmax,Imax); D.etaN= zeros(Jmax,Imax); D.etaSi = zeros(Jmax,Imax); 
G.muN = zeros(Jmax,Imax);D.muN = zeros(Jmax,Imax);D.muSi = zeros(Jmax,Imax); 

help1 = max(0,(1-betaL)*D.JL-D.Jresp-betaN*D.JNreg-betaSi*D.JSireg);
for i=1:Jmax*Imax
    if D.limind(i)==1 %C limiting
        D.etaN(i) = D.JNreg(i)-help1(i)./D.CN(i);
        D.etaSi(i) = D.JSireg(i)-help1(i)./D.CS(i);
    elseif D.limind(i)==2 %N limiting
        D.etaSi(i) = D.JSireg(i)-D.JNreg(i).*D.CN(i)./D.CS(i);
            %DetaN=0 par initialisation
    elseif D.limind(i)==3 %Si limiting
        D.etaN(i)=D.JNreg(i)-D.JSireg(i).*D.CS(i)./D.CN(i);
            %DetaSi=0 par initialisation
    elseif G.Jeff(i)==G.Ccons(i)
        G.etaN(i) = (G.Ncons(i)-G.Ccons(i))./G.CN(i); 
    elseif G.Jeff(i) ==0
        G.etaN(i) = G.JN(i); 
    end
end   
    G.etaN(G.etaN<0)=0; D.etaN(D.etaN<0)=0; D.etaSi(D.etaSi<0)=0; %N and Si in excess
    G.muN = (G.JN-G.etaN).*(ones(Jmax,Imax)-G.g).*Bpred; G.muN(G.muN<0)=0; 
    D.muN = (D.JN-D.etaN).*(ones(Jmax,Imax)-D.g).*Bpred; D.muN(D.muN<0)=0;
    D.muSi = (D.JSi-D.etaSi).*(ones(Jmax,Imax)-D.g).*Bpred; D.muSi(D.muSi<0)=0;
    
for i = 1:Imax  % calculated once
    D.m = D.m + alphaF*Bpred(i).*T(:,:,i).*(ones(Jmax,Imax)-D.g(:,i)); %predation by larger cells for diatoms 
 end

  G.m = D.m(1,:);
  D.m = palat*D.m; % mortality rates [day-1] due to internal predation 
  D.m(D.m<0)=0;G.m(G.m<0)=0;
  


%==== Differential equations ====%
dBdt_d = dilday.*D0 + (D.g.*Bd-D.mp.*Bd-D.m.*Bd -m2.*Bd.^2 -dilday.*Bd); %dynamic equation for diatoms' biomass %matrix Jmax x Imax
dBdt_g = dilday.*P0 + (G.g.*Bg-G.mp.*Bg-G.m.*Bg -m2.*Bg.^2 -dilday.*Bg);%dynamic equation for generalists' biomass
dBdt = [reshape(dBdt_d,[Jmax*Imax 1]);reshape(dBdt_g,[Jmax*Imax 1])];


if seasonal==1
    auxN = reminN*((m2*Bd+D.m).*Bd./D.CN + (m2*Bg+G.m).*Bg./G.CN) - D.muN - G.muN;
    auxSi = reminS*(m2*Bd+D.m).*Bd./D.CS - D.muSi;
    dBdt_N = dilday.*(N0-N) + sum(auxN(:)); %dynamic equation for nitrogen
    dBdt_Si = dilday.*(S0-Si) + sum(auxSi(:)); %dynamic equation for silica
    dBdt = [reshape(dBdt_d,[Jmax*Imax 1]);reshape(dBdt_g,[Jmax*Imax 1]); dBdt_N; dBdt_Si]; 
end
 
struct1=D; %output of structures D and G
struct2=G;
 
end
