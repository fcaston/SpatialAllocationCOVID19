clear all; close all

%Note that SIP (shelter-in-place) refers to travel restrictions.

%5 percent Drug
if true
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   Wrong Assumption About Immunity   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   5 Percent Drug   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% We assume immunity is temporary, while in fact it is permanent
%%% We assume correctly that there is compliance to SIP order
if true


load WorkspacePS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts05_D;

S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp'); %pchipinterp

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end


    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_SIP(x, beta11, beta22, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);

CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
DeltaCumI=(CumI_True-CumI_Optimal)/CumI_Optimal;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;


%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempInc_SIPCor_Drug05=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_TempInc_SIPCor_Drug05=DeltaNPV;
DeltaPenalty_TempInc_SIPCor_Drug05=DeltaPenalty;

end

%%% We assume immunity is temporary, while in fact it is permanent
%%% We assume correctly that there is noncompliance to SIP order
if true


load WorkspacePN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts05_D;

S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal));%/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempInc_noSIPCor_Drug05=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_TempInc_noSIPCor_Drug05=DeltaNPV;
DeltaPenalty_TempInc_noSIPCor_Drug05=DeltaPenalty;

end

%%% We assume immunity is permanent, while in fact it is temporary
%%% We assume correctly that there is compliance to SIP order
if true


load WorkspaceSS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts05_D;

S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspacePS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_SIP(x, beta11, beta22, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_PermInc_SIPCor_Drug05=DeltaCumI; 
DeltaNPV_PermInc_SIPCor_Drug05=DeltaNPV;
DeltaPenalty_PermInc_SIPCor_Drug05=DeltaPenalty;

end

%%% We assume immunity is permanent, while in fact it is temporary
%%% We assume correctly that there is noncompliance to SIP order
if true


load WorkspaceSN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts05_D;

S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspacePN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end


    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_PermInc_noSIPCor_Drug05=DeltaCumI; 
DeltaNPV_PermInc_noIPCor_Drug05=DeltaNPV;
DeltaPenalty_PermInc_noIPCor_Drug05=DeltaPenalty;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   Wrong Assumption About SIP   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   5 Percent Drug   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% We assume correctly that immunity is permanent
%%% We assume there is noncompliance to SIP order, while in fact there is compliance
if true


load WorkspacePS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts05_D;

S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp'); %pchipinterp

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspacePN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_SIP(x, beta11, beta22, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;


DeltaCumI_PermCor_noSIPInc_Drug05=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_PermCor_noSIPInc_Drug05=DeltaNPV;
DeltaPenalty_PermCor_noSIPInc_Drug05=DeltaPenalty;

end

%%% We assume correctly that immunity is permanent
%%% We assume there is compliance to SIP order, while in fact there is noncompliance
if true


load WorkspacePN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts05_D;

S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspacePS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;


DeltaCumI_PermCor_SIPInc_Drug05=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_PermCor_SIPInc_Drug05=DeltaNPV;
DeltaPenalty_PermCor_SIPInc_Drug05=DeltaPenalty;


end

%%% We assume correctly that immunity is temporary
%%% We assume there is noncompliance to SIP order, while in fact there is compliance
if true


load WorkspaceSS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts05_D;

S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end




    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_SIP(x, beta11, beta22, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempCor_noSIPInc_Drug05=DeltaCumI; 
DeltaNPV_TempCor_noSIPInc_Drug05=DeltaNPV;
DeltaPenalty_TempCor_noSIPInc_Drug05=DeltaPenalty;

end

%%% We assume correctly that immunity is temporary
%%% We assume there is compliance to SIP order, while in fact there is noncompliance
if true


load WorkspaceSN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts05_D;

S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*( cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;



DeltaCumI_TempCor_SIPInc_Drug05=DeltaCumI; 
DeltaNPV_TempCor_SIPInc_Drug05=DeltaNPV;
DeltaPenalty_TempCor_SIPInc_Drug05=DeltaPenalty;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Wrong Assumptions About Immunity and SIP   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   5 Percent Drug   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% We assume immunity is temporary, while in fact it is permanent
%%% We assume there is noncompliance to SIP order, while in fact there is compliance
if true


load WorkspacePS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts05_D;

S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp'); %pchipinterp

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_SIP(x, beta11, beta22, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempInc_noSIPInc_Drug05=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_TempInc_noSIPInc_Drug05=DeltaNPV;
DeltaPenalty_TempInc_noSIPInc_Drug05=DeltaPenalty;

end

%%% We assume immunity is temporary, while in fact it is permanent
%%% We assume there is compliance to SIP order, while in fact there is noncompliance
if true


load WorkspacePN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts05_D;

S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;




%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;



DeltaCumI_TempInc_SIPInc_Drug05=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_TempInc_SIPInc_Drug05=DeltaNPV;
DeltaPenalty_TempInc_SIPInc_Drug05=DeltaPenalty;

end

%%% We assume immunity is permanent, while in fact it is temporary
%%% We assume there is noncompliance to SIP order, while in fact there is compliance
if true


load WorkspaceSS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts05_D;

S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspacePN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end


    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_SIP(x, beta11, beta22, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;




%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;



DeltaCumI_PermInc_noSIPInc_Drug05=DeltaCumI; 
DeltaNPV_PermInc_noSIPInc_Drug05=DeltaNPV;
DeltaPenalty_PermInc_noSIPInc_Drug05=DeltaPenalty;

end

%%% We assume immunity is permanent, while in fact it is temporary
%%% We assume there is compliance to SIP order, while in fact there is noncompliance
if true


load WorkspaceSN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts05_D;

S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspacePS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;




%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;


%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;


DeltaCumI_PermInc_SIPInc_Drug05=DeltaCumI; 
DeltaNPV_PermInc_SIPInc_Drug05=DeltaNPV;
DeltaPenalty_PermInc_SIPInc_Drug05=DeltaPenalty;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Correct Assumptions but Ad Hoc   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   5 Percent Drug   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% We assume correctly that immunity is permanent
%%% We assume correctly that there is compliance to SIP order
%%% But we used the ad hoc allocation
if true


load WorkspacePS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts05_D;

S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp'); %pchipinterp

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);



end

%Comment below to see that the Deltas are equal to zero
load WorkspacePS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s05_Dah;
E1=E1s05_Dah;    
I1=I1s05_Dah;    
R1=R1s05_Dah;   
N1=N1s05_Dah;

S2=S2s05_Dah;
E2=E2s05_Dah;    
I2=I2s05_Dah;    
R2=R2s05_Dah;
N2=N2s05_Dah;

UD1=uD1s05_ah;
UD2=uD2s05_ah;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end


    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_SIP(x, beta11, beta22, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;




%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;


DeltaCumI_PermCor_SIPCor_Drug05AH=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_PermCor_SIPCor_Drug05AH=DeltaNPV;
DeltaPenalty_PermCor_SIPCor_Drug05AH=DeltaPenalty;

end


%%% We assume correctly that immunity is permanent
%%% We assume correctly that there is noncompliance to SIP order
%%% But we used the ad hoc allocation
if true


load WorkspacePN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts05_D;

S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspacePN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s05_Dah;
E1=E1s05_Dah;    
I1=I1s05_Dah;    
R1=R1s05_Dah;   
N1=N1s05_Dah;

S2=S2s05_Dah;
E2=E2s05_Dah;    
I2=I2s05_Dah;    
R2=R2s05_Dah;
N2=N2s05_Dah;

UD1=uD1s05_ah;
UD2=uD2s05_ah;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*( cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_PermCor_noSIPCor_Drug05AH=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_PermCor_noSIPCor_Drug05AH=DeltaNPV;
DeltaPenalty_PermCor_noSIPCor_Drug05AH=DeltaPenalty;

end


%%% We assume correctly that immunity is temporary
%%% We assume correctly that there is compliance to SIP order
%%% But we used the ad hoc allocation
if true


load WorkspaceSS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts05_D;

S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s05_Dah;
E1=E1s05_Dah;    
I1=I1s05_Dah;    
R1=R1s05_Dah;   
N1=N1s05_Dah;

S2=S2s05_Dah;
E2=E2s05_Dah;    
I2=I2s05_Dah;    
R2=R2s05_Dah;
N2=N2s05_Dah;

UD1=uD1s05_ah;
UD2=uD2s05_ah;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end


    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_SIP(x, beta11, beta22, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempCor_SIPCor_Drug05AH=DeltaCumI; 
DeltaNPV_TempCot_SIPCor_Drug05AH=DeltaNPV;
DeltaPenalty_TempCot_SIPCor_Drug05AH=DeltaPenalty;

end

%%% We assume correctly that immunity is temporary
%%% We assume correctly that there is noncompliance to SIP order
%%% But we used the ad hoc allocation
if true

load WorkspaceSN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts05_D;

S1=S1s05_D;
E1=E1s05_D;    
I1=I1s05_D;    
R1=R1s05_D;   
N1=N1s05_D;

S2=S2s05_D;
E2=E2s05_D;    
I2=I2s05_D;    
R2=R2s05_D;
N2=N2s05_D;

UD1=uD1s05;
UD2=uD2s05;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);


%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s05_Dah;
E1=E1s05_Dah;    
I1=I1s05_Dah;    
R1=R1s05_Dah;   
N1=N1s05_Dah;

S2=S2s05_Dah;
E2=E2s05_Dah;    
I2=I2s05_Dah;    
R2=R2s05_Dah;
N2=N2s05_Dah;

UD1=uD1s05_ah;
UD2=uD2s05_ah;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;


%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempCor_noSIPCor_Drug05AH=DeltaCumI; 
DeltaNPV_TempCor_noIPCor_Drug05AH=DeltaNPV;
DeltaPenalty_TempCor_noIPCor_Drug05AH=DeltaPenalty;

end


    
    
end

%10 percent Drug
if true
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   Wrong Assumption About Immunity   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   10 Percent Drug  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% We assume immunity is temporary, while in fact it is permanent
%%% We assume correctly that there is compliance to SIP order
if true


load WorkspacePS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts10_D;

S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp'); %pchipinterp

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end


    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_SIP(x, beta11, beta22, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);

CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
DeltaCumI=(CumI_True-CumI_Optimal)/CumI_Optimal;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;


%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempInc_SIPCor_Drug10=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_TempInc_SIPCor_Drug10=DeltaNPV;
DeltaPenalty_TempInc_SIPCor_Drug10=DeltaPenalty;

end

%%% We assume immunity is temporary, while in fact it is permanent
%%% We assume correctly that there is noncompliance to SIP order
if true


load WorkspacePN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts10_D;

S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal));%/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempInc_noSIPCor_Drug10=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_TempInc_noSIPCor_Drug10=DeltaNPV;
DeltaPenalty_TempInc_noSIPCor_Drug10=DeltaPenalty;

end

%%% We assume immunity is permanent, while in fact it is temporary
%%% We assume correctly that there is compliance to SIP order
if true


load WorkspaceSS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts10_D;

S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspacePS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_SIP(x, beta11, beta22, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_PermInc_SIPCor_Drug10=DeltaCumI; 
DeltaNPV_PermInc_SIPCor_Drug10=DeltaNPV;
DeltaPenalty_PermInc_SIPCor_Drug10=DeltaPenalty;

end

%%% We assume immunity is permanent, while in fact it is temporary
%%% We assume correctly that there is noncompliance to SIP order
if true


load WorkspaceSN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts10_D;

S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspacePN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end


    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_PermInc_noSIPCor_Drug10=DeltaCumI; 
DeltaNPV_PermInc_noIPCor_Drug10=DeltaNPV;
DeltaPenalty_PermInc_noIPCor_Drug10=DeltaPenalty;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   Wrong Assumption About SIP   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   10 Percent Drug  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% We assume correctly that immunity is permanent
%%% We assume there is noncompliance to SIP order, while in fact there is compliance
if true


load WorkspacePS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts10_D;

S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp'); %pchipinterp

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspacePN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_SIP(x, beta11, beta22, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;


DeltaCumI_PermCor_noSIPInc_Drug10=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_PermCor_noSIPInc_Drug10=DeltaNPV;
DeltaPenalty_PermCor_noSIPInc_Drug10=DeltaPenalty;

end

%%% We assume correctly that immunity is permanent
%%% We assume there is compliance to SIP order, while in fact there is noncompliance
if true


load WorkspacePN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts10_D;

S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspacePS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;


DeltaCumI_PermCor_SIPInc_Drug10=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_PermCor_SIPInc_Drug10=DeltaNPV;
DeltaPenalty_PermCor_SIPInc_Drug10=DeltaPenalty;


end

%%% We assume correctly that immunity is temporary
%%% We assume there is noncompliance to SIP order, while in fact there is compliance
if true


load WorkspaceSS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts10_D;

S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end




    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_SIP(x, beta11, beta22, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempCor_noSIPInc_Drug10=DeltaCumI; 
DeltaNPV_TempCor_noSIPInc_Drug10=DeltaNPV;
DeltaPenalty_TempCor_noSIPInc_Drug10=DeltaPenalty;

end

%%% We assume correctly that immunity is temporary
%%% We assume there is compliance to SIP order, while in fact there is noncompliance
if true


load WorkspaceSN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts10_D;

S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*( cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;



DeltaCumI_TempCor_SIPInc_Drug10=DeltaCumI; 
DeltaNPV_TempCor_SIPInc_Drug10=DeltaNPV;
DeltaPenalty_TempCor_SIPInc_Drug10=DeltaPenalty;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Wrong Assumptions About Immunity and SIP   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   10 Percent Drug  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% We assume immunity is temporary, while in fact it is permanent
%%% We assume there is noncompliance to SIP order, while in fact there is compliance
if true


load WorkspacePS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts10_D;

S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp'); %pchipinterp

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_SIP(x, beta11, beta22, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempInc_noSIPInc_Drug10=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_TempInc_noSIPInc_Drug10=DeltaNPV;
DeltaPenalty_TempInc_noSIPInc_Drug10=DeltaPenalty;

end

%%% We assume immunity is temporary, while in fact it is permanent
%%% We assume there is compliance to SIP order, while in fact there is noncompliance
if true


load WorkspacePN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts10_D;

S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;




%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;



DeltaCumI_TempInc_SIPInc_Drug10=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_TempInc_SIPInc_Drug10=DeltaNPV;
DeltaPenalty_TempInc_SIPInc_Drug10=DeltaPenalty;

end

%%% We assume immunity is permanent, while in fact it is temporary
%%% We assume there is noncompliance to SIP order, while in fact there is compliance
if true


load WorkspaceSS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts10_D;

S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspacePN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end


    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_SIP(x, beta11, beta22, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;




%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;



DeltaCumI_PermInc_noSIPInc_Drug10=DeltaCumI; 
DeltaNPV_PermInc_noSIPInc_Drug10=DeltaNPV;
DeltaPenalty_PermInc_noSIPInc_Drug10=DeltaPenalty;

end

%%% We assume immunity is permanent, while in fact it is temporary
%%% We assume there is compliance to SIP order, while in fact there is noncompliance
if true


load WorkspaceSN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts10_D;

S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspacePS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;




%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;


%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;


DeltaCumI_PermInc_SIPInc_Drug10=DeltaCumI; 
DeltaNPV_PermInc_SIPInc_Drug10=DeltaNPV;
DeltaPenalty_PermInc_SIPInc_Drug10=DeltaPenalty;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Correct Assumptions but Ad Hoc   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   10 Percent Drug  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% We assume correctly that immunity is permanent
%%% We assume correctly that there is compliance to SIP order
%%% But we used the ad hoc allocation
if true


load WorkspacePS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts10_D;

S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp'); %pchipinterp

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);



end

%Comment below to see that the Deltas are equal to zero
load WorkspacePS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s10_Dah;
E1=E1s10_Dah;    
I1=I1s10_Dah;    
R1=R1s10_Dah;   
N1=N1s10_Dah;

S2=S2s10_Dah;
E2=E2s10_Dah;    
I2=I2s10_Dah;    
R2=R2s10_Dah;
N2=N2s10_Dah;

UD1=uD1s10_ah;
UD2=uD2s10_ah;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end


    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_SIP(x, beta11, beta22, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;




%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;


DeltaCumI_PermCor_SIPCor_Drug10AH=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_PermCor_SIPCor_Drug10AH=DeltaNPV;
DeltaPenalty_PermCor_SIPCor_Drug10AH=DeltaPenalty;

end


%%% We assume correctly that immunity is permanent
%%% We assume correctly that there is noncompliance to SIP order
%%% But we used the ad hoc allocation
if true


load WorkspacePN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts10_D;

S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspacePN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s10_Dah;
E1=E1s10_Dah;    
I1=I1s10_Dah;    
R1=R1s10_Dah;   
N1=N1s10_Dah;

S2=S2s10_Dah;
E2=E2s10_Dah;    
I2=I2s10_Dah;    
R2=R2s10_Dah;
N2=N2s10_Dah;

UD1=uD1s10_ah;
UD2=uD2s10_ah;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*( cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_PermCor_noSIPCor_Drug10AH=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_PermCor_noSIPCor_Drug10AH=DeltaNPV;
DeltaPenalty_PermCor_noSIPCor_Drug10AH=DeltaPenalty;

end


%%% We assume correctly that immunity is temporary
%%% We assume correctly that there is compliance to SIP order
%%% But we used the ad hoc allocation
if true


load WorkspaceSS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts10_D;

S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s10_Dah;
E1=E1s10_Dah;    
I1=I1s10_Dah;    
R1=R1s10_Dah;   
N1=N1s10_Dah;

S2=S2s10_Dah;
E2=E2s10_Dah;    
I2=I2s10_Dah;    
R2=R2s10_Dah;
N2=N2s10_Dah;

UD1=uD1s10_ah;
UD2=uD2s10_ah;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end


    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_SIP(x, beta11, beta22, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempCor_SIPCor_Drug10AH=DeltaCumI; 
DeltaNPV_TempCot_SIPCor_Drug10AH=DeltaNPV;
DeltaPenalty_TempCot_SIPCor_Drug10AH=DeltaPenalty;

end

%%% We assume correctly that immunity is temporary
%%% We assume correctly that there is noncompliance to SIP order
%%% But we used the ad hoc allocation
if true

load WorkspaceSN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts10_D;

S1=S1s10_D;
E1=E1s10_D;    
I1=I1s10_D;    
R1=R1s10_D;   
N1=N1s10_D;

S2=S2s10_D;
E2=E2s10_D;    
I2=I2s10_D;    
R2=R2s10_D;
N2=N2s10_D;

UD1=uD1s10;
UD2=uD2s10;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);


%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s10_Dah;
E1=E1s10_Dah;    
I1=I1s10_Dah;    
R1=R1s10_Dah;   
N1=N1s10_Dah;

S2=S2s10_Dah;
E2=E2s10_Dah;    
I2=I2s10_Dah;    
R2=R2s10_Dah;
N2=N2s10_Dah;

UD1=uD1s10_ah;
UD2=uD2s10_ah;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;


%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempCor_noSIPCor_Drug10AH=DeltaCumI; 
DeltaNPV_TempCor_noIPCor_Drug10AH=DeltaNPV;
DeltaPenalty_TempCor_noIPCor_Drug10AH=DeltaPenalty;

end


    
    
end

%15 percent Drug
if true
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   Wrong Assumption About Immunity   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   15 Percent Drug  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% We assume immunity is temporary, while in fact it is permanent
%%% We assume correctly that there is compliance to SIP order
if true


load WorkspacePS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts15_D;

S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp'); %pchipinterp

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end


    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_SIP(x, beta11, beta22, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);

CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
DeltaCumI=(CumI_True-CumI_Optimal)/CumI_Optimal;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;


%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempInc_SIPCor_Drug15=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_TempInc_SIPCor_Drug15=DeltaNPV;
DeltaPenalty_TempInc_SIPCor_Drug15=DeltaPenalty;

end

%%% We assume immunity is temporary, while in fact it is permanent
%%% We assume correctly that there is noncompliance to SIP order
if true


load WorkspacePN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts15_D;

S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal));%/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempInc_noSIPCor_Drug15=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_TempInc_noSIPCor_Drug15=DeltaNPV;
DeltaPenalty_TempInc_noSIPCor_Drug15=DeltaPenalty;

end

%%% We assume immunity is permanent, while in fact it is temporary
%%% We assume correctly that there is compliance to SIP order
if true


load WorkspaceSS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts15_D;

S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspacePS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_SIP(x, beta11, beta22, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_PermInc_SIPCor_Drug15=DeltaCumI; 
DeltaNPV_PermInc_SIPCor_Drug15=DeltaNPV;
DeltaPenalty_PermInc_SIPCor_Drug15=DeltaPenalty;

end

%%% We assume immunity is permanent, while in fact it is temporary
%%% We assume correctly that there is noncompliance to SIP order
if true


load WorkspaceSN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts15_D;

S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspacePN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end


    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_PermInc_noSIPCor_Drug15=DeltaCumI; 
DeltaNPV_PermInc_noIPCor_Drug15=DeltaNPV;
DeltaPenalty_PermInc_noIPCor_Drug15=DeltaPenalty;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   Wrong Assumption About SIP   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   15 Percent Drug  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% We assume correctly that immunity is permanent
%%% We assume there is noncompliance to SIP order, while in fact there is compliance
if true


load WorkspacePS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts15_D;

S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp'); %pchipinterp

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspacePN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_SIP(x, beta11, beta22, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;


DeltaCumI_PermCor_noSIPInc_Drug15=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_PermCor_noSIPInc_Drug15=DeltaNPV;
DeltaPenalty_PermCor_noSIPInc_Drug15=DeltaPenalty;

end

%%% We assume correctly that immunity is permanent
%%% We assume there is compliance to SIP order, while in fact there is noncompliance
if true


load WorkspacePN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts15_D;

S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspacePS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;


DeltaCumI_PermCor_SIPInc_Drug15=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_PermCor_SIPInc_Drug15=DeltaNPV;
DeltaPenalty_PermCor_SIPInc_Drug15=DeltaPenalty;


end

%%% We assume correctly that immunity is temporary
%%% We assume there is noncompliance to SIP order, while in fact there is compliance
if true


load WorkspaceSS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts15_D;

S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end




    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_SIP(x, beta11, beta22, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempCor_noSIPInc_Drug15=DeltaCumI; 
DeltaNPV_TempCor_noSIPInc_Drug15=DeltaNPV;
DeltaPenalty_TempCor_noSIPInc_Drug15=DeltaPenalty;

end

%%% We assume correctly that immunity is temporary
%%% We assume there is compliance to SIP order, while in fact there is noncompliance
if true


load WorkspaceSN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts15_D;

S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*( cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;



DeltaCumI_TempCor_SIPInc_Drug15=DeltaCumI; 
DeltaNPV_TempCor_SIPInc_Drug15=DeltaNPV;
DeltaPenalty_TempCor_SIPInc_Drug15=DeltaPenalty;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Wrong Assumptions About Immunity and SIP   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   15 Percent Drug  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% We assume immunity is temporary, while in fact it is permanent
%%% We assume there is noncompliance to SIP order, while in fact there is compliance
if true


load WorkspacePS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts15_D;

S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp'); %pchipinterp

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_SIP(x, beta11, beta22, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempInc_noSIPInc_Drug15=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_TempInc_noSIPInc_Drug15=DeltaNPV;
DeltaPenalty_TempInc_noSIPInc_Drug15=DeltaPenalty;

end

%%% We assume immunity is temporary, while in fact it is permanent
%%% We assume there is compliance to SIP order, while in fact there is noncompliance
if true


load WorkspacePN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts15_D;

S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;




%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;



DeltaCumI_TempInc_SIPInc_Drug15=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_TempInc_SIPInc_Drug15=DeltaNPV;
DeltaPenalty_TempInc_SIPInc_Drug15=DeltaPenalty;

end

%%% We assume immunity is permanent, while in fact it is temporary
%%% We assume there is noncompliance to SIP order, while in fact there is compliance
if true


load WorkspaceSS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts15_D;

S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspacePN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end


    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_SIP(x, beta11, beta22, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;




%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;



DeltaCumI_PermInc_noSIPInc_Drug15=DeltaCumI; 
DeltaNPV_PermInc_noSIPInc_Drug15=DeltaNPV;
DeltaPenalty_PermInc_noSIPInc_Drug15=DeltaPenalty;

end

%%% We assume immunity is permanent, while in fact it is temporary
%%% We assume there is compliance to SIP order, while in fact there is noncompliance
if true


load WorkspaceSN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts15_D;

S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspacePS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;




%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;


%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;


DeltaCumI_PermInc_SIPInc_Drug15=DeltaCumI; 
DeltaNPV_PermInc_SIPInc_Drug15=DeltaNPV;
DeltaPenalty_PermInc_SIPInc_Drug15=DeltaPenalty;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Correct Assumptions but Ad Hoc   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   15 Percent Drug  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% We assume correctly that immunity is permanent
%%% We assume correctly that there is compliance to SIP order
%%% But we used the ad hoc allocation
if true


load WorkspacePS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts15_D;

S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp'); %pchipinterp

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);



end

%Comment below to see that the Deltas are equal to zero
load WorkspacePS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s15_Dah;
E1=E1s15_Dah;    
I1=I1s15_Dah;    
R1=R1s15_Dah;   
N1=N1s15_Dah;

S2=S2s15_Dah;
E2=E2s15_Dah;    
I2=I2s15_Dah;    
R2=R2s15_Dah;
N2=N2s15_Dah;

UD1=uD1s15_ah;
UD2=uD2s15_ah;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end


    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_SIP(x, beta11, beta22, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;




%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;


DeltaCumI_PermCor_SIPCor_Drug15AH=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_PermCor_SIPCor_Drug15AH=DeltaNPV;
DeltaPenalty_PermCor_SIPCor_Drug15AH=DeltaPenalty;

end


%%% We assume correctly that immunity is permanent
%%% We assume correctly that there is noncompliance to SIP order
%%% But we used the ad hoc allocation
if true


load WorkspacePN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts15_D;

S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspacePN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s15_Dah;
E1=E1s15_Dah;    
I1=I1s15_Dah;    
R1=R1s15_Dah;   
N1=N1s15_Dah;

S2=S2s15_Dah;
E2=E2s15_Dah;    
I2=I2s15_Dah;    
R2=R2s15_Dah;
N2=N2s15_Dah;

UD1=uD1s15_ah;
UD2=uD2s15_ah;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Permanent_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*( cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_PermCor_noSIPCor_Drug15AH=DeltaCumI; %Assumed Imm is Not Perm, while it is Perm, SIP
DeltaNPV_PermCor_noSIPCor_Drug15AH=DeltaNPV;
DeltaPenalty_PermCor_noSIPCor_Drug15AH=DeltaPenalty;

end


%%% We assume correctly that immunity is temporary
%%% We assume correctly that there is compliance to SIP order
%%% But we used the ad hoc allocation
if true


load WorkspaceSS_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts15_D;

S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);

%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);


end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSS_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s15_Dah;
E1=E1s15_Dah;    
I1=I1s15_Dah;    
R1=R1s15_Dah;   
N1=N1s15_Dah;

S2=S2s15_Dah;
E2=E2s15_Dah;    
I2=I2s15_Dah;    
R2=R2s15_Dah;
N2=N2s15_Dah;

UD1=uD1s15_ah;
UD2=uD2s15_ah;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end


    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_SIP(x, beta11, beta22, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;

%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempCor_SIPCor_Drug15AH=DeltaCumI; 
DeltaNPV_TempCot_SIPCor_Drug15AH=DeltaNPV;
DeltaPenalty_TempCot_SIPCor_Drug15AH=DeltaPenalty;

end

%%% We assume correctly that immunity is temporary
%%% We assume correctly that there is noncompliance to SIP order
%%% But we used the ad hoc allocation
if true

load WorkspaceSN_Gauss.mat %This is the true scenario
if true %Optimal 10% Drug 
ts=ts15_D;

S1=S1s15_D;
E1=E1s15_D;    
I1=I1s15_D;    
R1=R1s15_D;   
N1=N1s15_D;

S2=S2s15_D;
E2=E2s15_D;    
I2=I2s15_D;    
R2=R2s15_D;
N2=N2s15_D;

UD1=uD1s15;
UD2=uD2s15;
end 

%Initial Conditions
x0(:,1)=[S1(1),E1(1),I1(1),R1(1),S2(1),E2(1),I2(1),R2(1), N1(1), N2(1)]';
%Terminal conditions if we would have correctly guessed
Optimal_end = [S1(end),E1(end),I1(end),R1(end),S2(end),E2(end),I2(end),R2(end),N1(end),N2(end)]';
%Interpolation of data
if true
%Type of interpolation
fitApprox = fittype('pchipinterp');

%Interpolation of Cumulative cases
gI1_Optimal=fit(ts,(I1),fitApprox);
gI2_Optimal=fit(ts,(I2),fitApprox);

CumI1_Optimal=integrate(gI1_Optimal,T,0);
CumI2_Optimal=integrate(gI2_Optimal,T,0);

%Interpolation of NPV
gNPV_Optimal = fit(ts,(exp(-r.*ts).*((phi+w).*cI.*(I1+I2) + cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
NPV_Optimal=integrate(gNPV_Optimal,T,0);


%Interpolation of Penalty
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(UD1+UD2) + cAdj.*((I2./(I1+I2)).*UD1 - (I1./(I1+I2)).*UD2).^2)),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);

end

%Comment below to see that the Deltas are equal to zero
load WorkspaceSN_Gauss.mat % This is the assumed scenario
if true %Optimal 10% Drug, renaming everything assuming wrong scenario
S1=S1s15_Dah;
E1=E1s15_Dah;    
I1=I1s15_Dah;    
R1=R1s15_Dah;   
N1=N1s15_Dah;

S2=S2s15_Dah;
E2=E2s15_Dah;    
I2=I2s15_Dah;    
R2=R2s15_Dah;
N2=N2s15_Dah;

UD1=uD1s15_ah;
UD2=uD2s15_ah;
end 

%Time steps
ts=ts10; %Time steps used to find the optimal solution
t_rk=0:(1/10000):T; %Runge-Kutta time steps

%Interpolation of control variables
uD1=interp1(ts,UD1,t_rk,'pchip');
uD2=interp1(ts,UD2,t_rk,'pchip');

%Simulating the optimal solution of the assumed scenario with the correct ODEs
omega=1/6;
for i = 2:length(t_rk+1)


if uD1(i-1) > x0(3,i-1)
    
uD1(i-1)=x0(3,i-1);

else

uD1(i-1)=uD1(i-1);

end

if uD2(i-1) > x0(7,i-1)
    
uD2(i-1)=x0(7,i-1);

else

uD2(i-1)=uD2(i-1);

end



    [t,y] = ode45(@(t,x)COVIDeqs_Temporary_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, omega, phi, qV, uD1(i-1), uD2(i-1)),[t_rk(i-1) t_rk(i)],x0(:,i-1));
    x0(:,i) = y(end,:);                 
end

%Storing the results
Results_RK=[x0;uD1;uD2];
%Naming the state variables that we obtained given we incorrectly guessed
if true
   
    S1_True=Results_RK(1,:)';
    E1_True=Results_RK(2,:)';
    I1_True=Results_RK(3,:)';
    R1_True=Results_RK(4,:)';
    S2_True=Results_RK(5,:)';
    E2_True=Results_RK(6,:)';
    I2_True=Results_RK(7,:)';
    R2_True=Results_RK(8,:)';
    N1_True=Results_RK(9,:)';
    N2_True=Results_RK(10,:)';
    
end

%Terminal conditions given we incorrectly guessed
True_end= x0(:,end);

%Difference in terminal conditions
Delta_FinalStates=[True_end-Optimal_end] % This should be approximately zero when correctly guess

%Interpolation of Cumulative cases resulting from the incorrect guess
gI1_True=fit(t_rk',(I1_True),fitApprox);
gI2_True=fit(t_rk',(I2_True),fitApprox);

%Calculating differences in infections
CumI1_True=integrate(gI1_True,T,0);
CumI2_True=integrate(gI2_True,T,0);
CumI_True= CumI1_True + CumI2_True; 
CumI_Optimal= CumI1_Optimal+CumI2_Optimal;
pop=1e0;
DeltaCumI=((CumI_True-CumI_Optimal)/CumI_Optimal).*pop;



%Interpolation of NPV
t_rk=t_rk'; uD1=uD1'; uD2=uD2';
gNPV_True = fit(t_rk,(exp(-r.*t_rk).*((phi+w).*cI.*(I1_True+I2_True) + cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
NPV_True=integrate(gNPV_True,T,0) ;
DeltaNPV=(NPV_True - NPV_Optimal)/NPV_Optimal;


%Interpolation of Penalty
gPenalty_True = fit(t_rk,(exp(-r.*t_rk).*(cD.*(uD1+uD2) + cAdj.*((I2_True./(I1_True+I2_True)).*uD1 - (I1_True./(I1_True+I2_True)).*uD2).^2)),fitApprox); 
Penalty_True=integrate(gPenalty_True,T,0) ;
DeltaPenalty=(Penalty_True - Penalty_Optimal)/Penalty_Optimal;

DeltaCumI_TempCor_noSIPCor_Drug15AH=DeltaCumI; 
DeltaNPV_TempCor_noIPCor_Drug15AH=DeltaNPV;
DeltaPenalty_TempCor_noIPCor_Drug15AH=DeltaPenalty;

end


    
    
end


%Figure with Penalty four panels 5%
if true


fig=figure
subplot(221)
plot(DeltaCumI_TempInc_SIPCor_Drug05*100,DeltaPenalty_TempInc_SIPCor_Drug05*100,'k+','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_PermCor_noSIPInc_Drug05*100,DeltaPenalty_PermCor_noSIPInc_Drug05*100,'ko','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_TempInc_noSIPInc_Drug05*100,DeltaPenalty_TempInc_noSIPInc_Drug05*100,'k*','LineWidth',0.75,'MarkerSize', 12); hold on
plot(DeltaCumI_PermCor_SIPCor_Drug05AH*100,DeltaPenalty_PermCor_SIPCor_Drug05AH*100,'k^','LineWidth',1,'MarkerSize', 12); hold on
ylim([-.35*100 .35*100])
xlim([0 0.04*100])
hline = refline([0 0]);
hline.Color = 'k';
hline.LineStyle = ':';
hline.HandleVisibility = 'off';
title({'Permanent Immunity &','Compliance to Travel Restrictions','(A)'})


subplot(222)
p1=plot(DeltaCumI_TempInc_noSIPCor_Drug05*100,DeltaPenalty_TempInc_noSIPCor_Drug05*100,'k+','LineWidth',1,'MarkerSize', 12); hold on
p2=plot(DeltaCumI_PermCor_SIPInc_Drug05*100,DeltaPenalty_PermCor_SIPInc_Drug05*100,'ko','LineWidth',1,'MarkerSize', 12); hold on
p3=plot(DeltaCumI_TempInc_SIPInc_Drug05*100,DeltaPenalty_TempInc_SIPInc_Drug05*100,'k*','LineWidth',0.75,'MarkerSize', 12); hold on
p4=plot(DeltaCumI_PermCor_noSIPCor_Drug05AH*100,DeltaPenalty_PermCor_noSIPCor_Drug05AH*100,'k^','LineWidth',1,'MarkerSize', 12); hold on
ylim([-.35*100 .35*100])
xlim([0 0.04*100])
hline = refline([0 0]);
hline.Color = 'k';
hline.LineStyle = ':';
hline.HandleVisibility = 'off';
title({'Permanent Immunity &','Noncompliance to Travel Restrictions','(B)'})
      legend1=legend([p1 p3 p2 p4],{'Immunity','Travel Restrictions','Both','Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.219215899371855 0.483039942998033 0.582857693944659 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

   


subplot(223)
plot(DeltaCumI_PermInc_SIPCor_Drug05*100,DeltaPenalty_PermInc_SIPCor_Drug05*100,'k+','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_TempCor_noSIPInc_Drug05*100,DeltaPenalty_TempCor_noSIPInc_Drug05*100,'ko','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_PermInc_noSIPInc_Drug05*100,DeltaPenalty_PermInc_noSIPInc_Drug05*100,'k*','LineWidth',0.75,'MarkerSize', 12); hold on
plot(DeltaCumI_TempCor_SIPCor_Drug05AH*100,DeltaPenalty_TempCot_SIPCor_Drug05AH*100,'k^','LineWidth',1,'MarkerSize', 12); hold on
ylim([-.35*100 .35*100])
xlim([0 0.04*100])
hline = refline([0 0]);
hline.Color = 'k';
hline.LineStyle = ':';
hline.HandleVisibility = 'off';
title({'','6-month Immunity &','Compliance to Travel Restrictions','(C)'})

subplot(224)
plot(DeltaCumI_PermInc_noSIPCor_Drug05*100,DeltaPenalty_PermInc_noIPCor_Drug05*100,'k+','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_TempCor_SIPInc_Drug05*100,DeltaPenalty_TempCor_SIPInc_Drug05*100,'ko','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_PermInc_SIPInc_Drug05*100,DeltaPenalty_PermInc_SIPInc_Drug05*100,'k*','LineWidth',0.75,'MarkerSize', 12); hold on
plot(DeltaCumI_TempCor_noSIPCor_Drug05AH*100,DeltaPenalty_TempCor_noIPCor_Drug05AH*100,'k^','LineWidth',1,'MarkerSize', 12); hold on
ylim([-.35*100 .35*100])
xlim([0 0.04*100])
title({'','6-month Immunity &','Noncompliance to Travel Restrictions','(D)'})
hline = refline([0 0]);
hline.Color = 'k';
hline.LineStyle = ':';
hline.HandleVisibility = 'off';


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 
    xlabel(han, {'Percentage Change in Cumulative Cases'}, 'FontSize', 16);
    ylabel(han, {'Percentage Change in Expenditures'}, 'FontSize', 16); 

%sgtitle({'Cross: wrong about immunity length; Circle: wrong about SIP','Star: wrong about both; Triangle: Ad Hoc','Blue: Perm. and SIP is true','Green: Perm and no SIP is true','Red: Temporary and SIP is true','Black: Temporary and no SIP is true'})



   saveas(gcf,'WrongAssumption_Drugs_Total_5.png'); hold off
end


%Figure with Penalty four panels 10%
if true


fig=figure
subplot(221)
plot(DeltaCumI_TempInc_SIPCor_Drug10*100,DeltaPenalty_TempInc_SIPCor_Drug10*100,'k+','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_PermCor_noSIPInc_Drug10*100,DeltaPenalty_PermCor_noSIPInc_Drug10*100,'ko','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_TempInc_noSIPInc_Drug10*100,DeltaPenalty_TempInc_noSIPInc_Drug10*100,'k*','LineWidth',0.75,'MarkerSize', 12); hold on
plot(DeltaCumI_PermCor_SIPCor_Drug10AH*100,DeltaPenalty_PermCor_SIPCor_Drug10AH*100,'k^','LineWidth',1,'MarkerSize', 12); hold on
ylim([-.35*100 .35*100])
xlim([0 0.04*100])
hline = refline([0 0]);
hline.Color = 'k';
hline.LineStyle = ':';
hline.HandleVisibility = 'off';
title({'Permanent Immunity &','Compliance to Travel Restrictions','(A)'})


subplot(222)
p1=plot(DeltaCumI_TempInc_noSIPCor_Drug10*100,DeltaPenalty_TempInc_noSIPCor_Drug10*100,'k+','LineWidth',1,'MarkerSize', 12); hold on
p2=plot(DeltaCumI_PermCor_SIPInc_Drug10*100,DeltaPenalty_PermCor_SIPInc_Drug10*100,'ko','LineWidth',1,'MarkerSize', 12); hold on
p3=plot(DeltaCumI_TempInc_SIPInc_Drug10*100,DeltaPenalty_TempInc_SIPInc_Drug10*100,'k*','LineWidth',0.75,'MarkerSize', 12); hold on
p4=plot(DeltaCumI_PermCor_noSIPCor_Drug10AH*100,DeltaPenalty_PermCor_noSIPCor_Drug10AH*100,'k^','LineWidth',1,'MarkerSize', 12); hold on
ylim([-.35*100 .35*100])
xlim([0 0.04*100])
hline = refline([0 0]);
hline.Color = 'k';
hline.LineStyle = ':';
hline.HandleVisibility = 'off';
title({'Permanent Immunity &','Noncompliance to Travel Restrictions','(B)'})
  legend1=legend([p1 p3 p2 p4],{'Immunity','Travel Restrictions','Both','Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.219215899371855 0.483039942998033 0.582857693944659 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

   


subplot(223)
plot(DeltaCumI_PermInc_SIPCor_Drug10*100,DeltaPenalty_PermInc_SIPCor_Drug10*100,'k+','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_TempCor_noSIPInc_Drug10*100,DeltaPenalty_TempCor_noSIPInc_Drug10*100,'ko','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_PermInc_noSIPInc_Drug10*100,DeltaPenalty_PermInc_noSIPInc_Drug10*100,'k*','LineWidth',0.75,'MarkerSize', 12); hold on
plot(DeltaCumI_TempCor_SIPCor_Drug10AH*100,DeltaPenalty_TempCot_SIPCor_Drug10AH*100,'k^','LineWidth',1,'MarkerSize', 12); hold on
ylim([-.35*100 .35*100])
xlim([0 0.04*100])
hline = refline([0 0]);
hline.Color = 'k';
hline.LineStyle = ':';
hline.HandleVisibility = 'off';
title({'','6-month Immunity &','Compliance to Travel Restrictions','(C)'})

subplot(224)
plot(DeltaCumI_PermInc_noSIPCor_Drug10*100,DeltaPenalty_PermInc_noIPCor_Drug10*100,'k+','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_TempCor_SIPInc_Drug10*100,DeltaPenalty_TempCor_SIPInc_Drug10*100,'ko','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_PermInc_SIPInc_Drug10*100,DeltaPenalty_PermInc_SIPInc_Drug10*100,'k*','LineWidth',0.75,'MarkerSize', 12); hold on
plot(DeltaCumI_TempCor_noSIPCor_Drug10AH*100,DeltaPenalty_TempCor_noIPCor_Drug10AH*100,'k^','LineWidth',1,'MarkerSize', 12); hold on
ylim([-.35*100 .35*100])
xlim([0 0.04*100])
title({'','6-month Immunity &','Noncompliance to Travel Restrictions','(D)'})
hline = refline([0 0]);
hline.Color = 'k';
hline.LineStyle = ':';
hline.HandleVisibility = 'off';


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 
    xlabel(han, {'Percentage Change in Cumulative Cases'}, 'FontSize', 16);
    ylabel(han, {'Percentage Change in Expenditures'}, 'FontSize', 16); 

%sgtitle({'Cross: wrong about immunity length; Circle: wrong about SIP','Star: wrong about both; Triangle: Ad Hoc','Blue: Perm. and SIP is true','Green: Perm and no SIP is true','Red: Temporary and SIP is true','Black: Temporary and no SIP is true'})



   saveas(gcf,'WrongAssumption_Drugs_Total_10.png'); hold off
end


%Figure with Penalty four panels 15%
if true


fig=figure
subplot(221)
plot(DeltaCumI_TempInc_SIPCor_Drug15*100,DeltaPenalty_TempInc_SIPCor_Drug15*100,'k+','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_PermCor_noSIPInc_Drug15*100,DeltaPenalty_PermCor_noSIPInc_Drug15*100,'ko','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_TempInc_noSIPInc_Drug15*100,DeltaPenalty_TempInc_noSIPInc_Drug15*100,'k*','LineWidth',0.75,'MarkerSize', 12); hold on
plot(DeltaCumI_PermCor_SIPCor_Drug15AH*100,DeltaPenalty_PermCor_SIPCor_Drug15AH*100,'k^','LineWidth',1,'MarkerSize', 12); hold on
ylim([-.45*100 .45*100])
xlim([0 0.11*100])
hline = refline([0 0]);
hline.Color = 'k';
hline.LineStyle = ':';
hline.HandleVisibility = 'off';
title({'Permanent Immunity &','Compliance to Travel Restrictions','(A)'})


subplot(222)
p1=plot(DeltaCumI_TempInc_noSIPCor_Drug15*100,DeltaPenalty_TempInc_noSIPCor_Drug15*100,'k+','LineWidth',1,'MarkerSize', 12); hold on
p2=plot(DeltaCumI_PermCor_SIPInc_Drug15*100,DeltaPenalty_PermCor_SIPInc_Drug15*100,'ko','LineWidth',1,'MarkerSize', 12); hold on
p3=plot(DeltaCumI_TempInc_SIPInc_Drug15*100,DeltaPenalty_TempInc_SIPInc_Drug15*100,'k*','LineWidth',0.75,'MarkerSize', 12); hold on
p4=plot(DeltaCumI_PermCor_noSIPCor_Drug15AH*100,DeltaPenalty_PermCor_noSIPCor_Drug15AH*100,'k^','LineWidth',1,'MarkerSize', 12); hold on
ylim([-.45*100 .45*100])
xlim([0 0.11*100])
hline = refline([0 0]);
hline.Color = 'k';
hline.LineStyle = ':';
hline.HandleVisibility = 'off';
title({'Permanent Immunity &','Noncompliance to Travel Restrictions','(B)'})
      legend1=legend([p1 p3 p2 p4],{'Immunity','Travel Restrictions','Both','Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
   'Position',[0.219215899371855 0.483039942998033 0.582857693944659 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

   


subplot(223)
plot(DeltaCumI_PermInc_SIPCor_Drug15*100,DeltaPenalty_PermInc_SIPCor_Drug15*100,'k+','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_TempCor_noSIPInc_Drug15*100,DeltaPenalty_TempCor_noSIPInc_Drug15*100,'ko','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_PermInc_noSIPInc_Drug15*100,DeltaPenalty_PermInc_noSIPInc_Drug15*100,'k*','LineWidth',0.75,'MarkerSize', 12); hold on
plot(DeltaCumI_TempCor_SIPCor_Drug15AH*100,DeltaPenalty_TempCot_SIPCor_Drug15AH*100,'k^','LineWidth',1,'MarkerSize', 12); hold on
ylim([-.45*100 .45*100])
xlim([0 0.11*100])
hline = refline([0 0]);
hline.Color = 'k';
hline.LineStyle = ':';
hline.HandleVisibility = 'off';
title({'','6-month Immunity &','Compliance to Travel Restrictions','(C)'})

subplot(224)
plot(DeltaCumI_PermInc_noSIPCor_Drug15*100,DeltaPenalty_PermInc_noIPCor_Drug15*100,'k+','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_TempCor_SIPInc_Drug15*100,DeltaPenalty_TempCor_SIPInc_Drug15*100,'ko','LineWidth',1,'MarkerSize', 12); hold on
plot(DeltaCumI_PermInc_SIPInc_Drug15*100,DeltaPenalty_PermInc_SIPInc_Drug15*100,'k*','LineWidth',0.75,'MarkerSize', 12); hold on
plot(DeltaCumI_TempCor_noSIPCor_Drug15AH*100,DeltaPenalty_TempCor_noIPCor_Drug15AH*100,'k^','LineWidth',1,'MarkerSize', 12); hold on
ylim([-.45*100 .45*100])
xlim([0 0.11*100])
title({'','6-month Immunity &','Noncompliance to Travel Restrictions','(D)'})
hline = refline([0 0]);
hline.Color = 'k';
hline.LineStyle = ':';
hline.HandleVisibility = 'off';


    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 
    xlabel(han, {'Percentage Change in Cumulative Cases'}, 'FontSize', 16);
    ylabel(han, {'Percentage Change in Expenditures'}, 'FontSize', 16); 

%sgtitle({'Cross: wrong about immunity length; Circle: wrong about SIP','Star: wrong about both; Triangle: Ad Hoc','Blue: Perm. and SIP is true','Green: Perm and no SIP is true','Red: Temporary and SIP is true','Black: Temporary and no SIP is true'})



   saveas(gcf,'WrongAssumption_Drugs_Total_15.png'); hold off
end




        function dx=COVIDeqs_Permanent_SIP(x, beta11, beta22, gamma, sigma, phi, qD, uD1, uD2)
        dx = zeros(10,1);

        dx(1) = - beta11.*x(1).*(x(3)./x(9)) ; 
        dx(2) = beta11.*x(1).*(x(3)./x(9)) - sigma.*x(2)  ; 
        dx(3) = sigma.*x(2) - gamma.*x(3) - phi.*x(3) - qD.*uD1 ; 
        dx(4) = gamma.*x(3)  + qD.*uD1;   
        dx(5) = - beta22.*x(5).*(x(7)./x(10))  ; 
        dx(6) = beta22.*x(5).*(x(7)./x(10)) - sigma.*x(6)  ; 
        dx(7) = sigma.*x(6) - gamma.*x(7) - phi.*x(7) - qD.*uD2 ; 
        dx(8) = gamma.*x(7)  + qD.*uD2;    
        dx(9) = - phi.*x(3);
        dx(10)= - phi.*x(7);
 
        end
    
        
        function dx=COVIDeqs_Permanent_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, phi, qD, uD1, uD2)
        dx = zeros(10,1);

        dx(1) = - beta11.*x(1).*(x(3)./x(9)) - beta12.*x(1).*(x(7)./x(10)) ; 
        dx(2) = beta11.*x(1).*(x(3)./x(9)) + beta12.*x(1).*(x(7)./x(10)) - sigma.*x(2)  ; 
        dx(3) = sigma.*x(2) - gamma.*x(3) - phi.*x(3) - qD.*uD1 ; 
        dx(4) = gamma.*x(3)  + qD.*uD1;   
        dx(5) = - beta22.*x(5).*(x(7)./x(10)) - beta21.*x(5).*(x(3)./x(9))  ; 
        dx(6) = beta22.*x(5).*(x(7)./x(10)) + beta21.*x(5).*(x(3)./x(9)) - sigma.*x(6)  ; 
        dx(7) = sigma.*x(6) - gamma.*x(7) - phi.*x(7) - qD.*uD2 ; 
        dx(8) = gamma.*x(7)  + qD.*uD2;    
        dx(9) = - phi.*x(3);
        dx(10)= - phi.*x(7);
 
        end
    
        
        function dx=COVIDeqs_Temporary_SIP(x, beta11, beta22, gamma, sigma, omega, phi, qD, uD1, uD2)
        dx = zeros(10,1);

        dx(1) = omega.*x(4) - beta11.*x(1).*(x(3)./x(9)) ; 
        dx(2) = beta11.*x(1).*(x(3)./x(9)) - sigma.*x(2)  ; 
        dx(3) = sigma.*x(2) - gamma.*x(3) - phi.*x(3) - qD.*uD1 ; 
        dx(4) = gamma.*x(3)  + qD.*uD1 - omega.*x(4);   
        dx(5) = omega.*x(8) - beta22.*x(5).*(x(7)./x(10))  ; 
        dx(6) = beta22.*x(5).*(x(7)./x(10)) - sigma.*x(6)  ; 
        dx(7) = sigma.*x(6) - gamma.*x(7) - phi.*x(7) - qD.*uD2 ; 
        dx(8) = gamma.*x(7)  + qD.*uD2 - omega.*x(8);    
        dx(9) = - phi.*x(3);
        dx(10)= - phi.*x(7);
 
        end
    
        
        function dx=COVIDeqs_Temporary_noSIP(x, beta11, beta22, beta12, beta21, gamma, sigma, omega, phi, qD, uD1, uD2)
        dx = zeros(10,1);

        dx(1) = omega.*x(4) - beta11.*x(1).*(x(3)./x(9)) - beta12.*x(1).*(x(7)./x(10)) ; 
        dx(2) = beta11.*x(1).*(x(3)./x(9)) + beta12.*x(1).*(x(7)./x(10)) - sigma.*x(2)  ; 
        dx(3) = sigma.*x(2) - gamma.*x(3) - phi.*x(3) - qD.*uD1 ; 
        dx(4) = gamma.*x(3)  + qD.*uD1 - omega.*x(4);   
        dx(5) = omega.*x(8) - beta22.*x(5).*(x(7)./x(10)) - beta21.*x(5).*(x(3)./x(9))  ; 
        dx(6) = beta22.*x(5).*(x(7)./x(10)) + beta21.*x(5).*(x(3)./x(9)) - sigma.*x(6)  ; 
        dx(7) = sigma.*x(6) - gamma.*x(7) - phi.*x(7) - qD.*uD2 ; 
        dx(8) = gamma.*x(7)  + qD.*uD2 - omega.*x(8);    
        dx(9) = - phi.*x(3);
        dx(10)= - phi.*x(7);
 
        end
    