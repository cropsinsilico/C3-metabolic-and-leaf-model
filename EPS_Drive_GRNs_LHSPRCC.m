function Arate=EPS_Drive_GRNs(input1,CO2i,PPFDi,Tempi,Gc,GT,Einput,k1,k2,k3,k4,c1,c2,c3,sens_para_vec)
global GRNC;
global GRNT;
global pcfactor;
global alpha1 alpha2;

GRN_data=input1*pcfactor;

CO2in = CO2i;
Liin = PPFDi;
Tpin = Tempi;
GRNC = Gc;
GRNT = GT;

global EnzymeAct;
EnzymeAct    = Einput(1:27)/30;%unit change
EnzymeAct(1) = EnzymeAct(1) * alpha1;
%This is part of scaling on the electron transport
%A few others are in CM_Rate.m
EnzymeAct(2:11) = EnzymeAct(2:11) * alpha2;
%Yufeng: test FBP/SBPase overexpression
EnzymeAct(5) = EnzymeAct(5) * 1; %V6:FBPase
EnzymeAct(8) = EnzymeAct(8) * 1; %V9:SBPase
%Yufeng: multiple a vector of scaling factors for all enzymes to be used in LHS-PRCC
% disp(size(sens_para_vec))
% disp(size(EnzymeAct))
EnzymeAct(1:26) = EnzymeAct(1:26) .* sens_para_vec'; %need a transpose here!
global Jmax;
Jmax=EnzymeAct(27);
global BFVmax;
BFVmax=Einput(28:45);
global FIVmax;
FIVmax=Einput(46:66);

global cATPsyn;
global CPSi;
global cNADPHsyn;
global cpsii;
if GRNC==1
cATPsyn=GRN_data(34);%1.0447;%1.01866 WY201803
CPSi=GRN_data(35);%1.0131;% 1.0237 WY201803
cNADPHsyn=GRN_data(37);%1.094468408;%1.0388 WY201803
cpsii=GRN_data(36);%1.0169;% 1.0129;%WY201803
end
if GRNC==0
cATPsyn=1;%1.01866 WY201803
CPSi=1;% 1.0237 WY201803
cNADPHsyn=1;%1.0388 WY201803
cpsii=1;% 1.0129;%WY201803   
end
global VfactorC;
global VfactorT;

VfactorC=GRN_data(1:33);

Arate=EPS_Drive(Liin,CO2in,Tpin,k1,k2,k3,k4,c1,c2,c3);
end
