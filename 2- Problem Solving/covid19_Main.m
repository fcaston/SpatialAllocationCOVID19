clear all; close all

%Calling parameter values
[N1,N2,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,r,cD,cV,cAdj,cI,w,params,R0] = covid_Parameters_months()

% Number of time periods
T =4; %(in months)

%Number of collocation points
Nset= 60; 
 
%Objective Function
OBJ=1; %=1 Minimize damages & costs

%Initial conditions 
x0ic=  [0.9074*N1; 0.0103*N1; 0.0143*N1; 0.0667*N1; 0.8662*N2; 0.0138*N2; 0.0196*N2; 0.0986*N2; N1; N2; 0; 0];

%%%Key varying factors
% The user needs to select which scenario
omega=0  ; ODE=1; %Permanent immunity, compliance to TR
%omega=0  ; ODE=2; %Permanent immunity, no compliance to TR
%omega=1/6; ODE=1; %6-month immunity, compliance to TR
%omega=1/6; ODE=2; %6-month immunity, no compliance to TR 

%Maximum phyical constraint (does not account for potential scarcity) 
MaxTreat=[0.1.*(N1+N2);N1+N2]; %constraint for uncontrolled R0 of 2.2 (it's 0.1091 in fact)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Main Results   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%No control case
if false
    
CASE=1;
    
SCAR=1;
    
%Initial guess to find unconstained amount
GUESS=[]; % = empty : Initial guess is initial conditions
    
omega=0; ODE=1; %Permanent immunity, compliance to TR 
    
    [Results, solution, ts, S1s_Dno1PS, S2s_Dno1PS, E1s_Dno1PS, E2s_Dno1PS, I1s_Dno1PS, I2s_Dno1PS, R1s_Dno1PS, R2s_Dno1PS, N1s_Dno1PS, N2s_Dno1PS, uD1s, uD2s, D1s, D2s] = ...
        covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

    [Results, solution, ts, S1s_no1PS, S2s_no1PS, E1s_no1PS, E2s_no1PS, I1s_no1PS, I2s_no1PS, R1s_no1PS, R2s_no1PS, N1s_no1PS, N2s_no1PS, uV1s, uV2s, V1s, V2s] = ...
        covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

omega=0; ODE=2;%6-month immunity, compliance to TR 
    
    [Results, solution, ts, S1s_Dno3PN, S2s_Dno3PN, E1s_Dno3PN, E2s_Dno3PN, I1s_Dno3PN, I2s_Dno3PN, R1s_Dno3PN, R2s_Dno3PN, N1s_Dno3PN, N2s_Dno3PN, uD1s, uD2s, D1s, D2s] = ...
        covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

    [Results, solution, ts, S1s_no3PN, S2s_no3PN, E1s_no3PN, E2s_no3PN, I1s_no3PN, I2s_no3PN, R1s_no3PN, R2s_no3PN, N1s_no3PN, N2s_no3PN, uV1s, uV2s, V1s, V2s] = ...
        covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);
    
omega=1/6; ODE=1; 
    
    [Results, solution, ts, S1s_Dno2SS, S2s_Dno2SS, E1s_Dno2SS, E2s_Dno2SS, I1s_Dno2SS, I2s_Dno2SS, R1s_Dno2SS, R2s_Dno2SS, N1s_Dno2SS, N2s_Dno2SS, uD1s, uD2s, D1s, D2s] = ...
        covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

    [Results, solution, ts, S1s_no2SS, S2s_no2SS, E1s_no2SS, E2s_no2SS, I1s_no2SS, I2s_no2SS, R1s_no2SS, R2s_no2SS, N1s_no2SS, N2s_no2SS, uV1s, uV2s, V1s, V2s] = ...
        covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);
    
omega=1/6; ODE=2; %6-month immunity, no compliance to TR
    
        [Results, solution, ts, S1s_Dno4SN, S2s_Dno4SN, E1s_Dno4SN, E2s_Dno4SN, I1s_Dno4SN, I2s_Dno4SN, R1s_Dno4SN, R2s_Dno4SN, N1s_Dno4SN, N2s_Dno4SN, uD1s, uD2s, D1s, D2s] = ...
        covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

    [Results, solution, ts, S1s_no4SN, S2s_no4SN, E1s_no4SN, E2s_no4SN, I1s_no4SN, I2s_no4SN, R1s_no4SN, R2s_no4SN, N1s_no4SN, N2s_no4SN, uV1s, uV2s, V1s, V2s] = ...
        covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);
    
    
end

%Optimal and Ad Hoc Allocations of Vaccines and Drugs
if false
    
    
    
%%% Allocation of antiviral drugs under different drug scarcity
if true

    %Initial guess to find unconstained amount
GUESS=[]; % = empty : Initial guess is initial conditions

CASE=2;



SCAR=0.025;

    
[Results, solution, ts05_D, S1s05_D, S2s05_D, E1s05_D, E2s05_D, I1s05_D, I2s05_D, R1s05_D, R2s05_D, N1s05_D, N2s05_D, uD1s05, uD2s05, D1s05, D2s05] = ...
 covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

if true %Uses results from SCAR=0.1 for a guess for SCAR=0.15
    GUESS(2).S1=S1s05_D;
    GUESS(2).E1=E1s05_D;
    GUESS(2).I1=I1s05_D;
    GUESS(2).R1=R1s05_D;
    GUESS(2).N1=N1s05_D;
    GUESS(2).S2=S2s05_D;
    GUESS(2).E2=E2s05_D;
    GUESS(2).I2=I2s05_D;
    GUESS(2).R2=R2s05_D;
    GUESS(2).N2=N2s05_D;
    GUESS(2).uD1=uD1s05*2;
    GUESS(2).uD2=uD2s05*2;
    GUESS(2).D1=D1s05;
    GUESS(2).D2=D2s05;
end

SCAR=0.05;
[Results, solution, ts05_D, S1s05_D, S2s05_D, E1s05_D, E2s05_D, I1s05_D, I2s05_D, R1s05_D, R2s05_D, N1s05_D, N2s05_D, uD1s05, uD2s05, D1s05, D2s05] = ...
    covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);



if true %Uses results from SCAR=0.1 for a guess for SCAR=0.15
    GUESS(2).S1=S1s05_D;
    GUESS(2).E1=E1s05_D;
    GUESS(2).I1=I1s05_D;
    GUESS(2).R1=R1s05_D;
    GUESS(2).N1=N1s05_D;
    GUESS(2).S2=S2s05_D;
    GUESS(2).E2=E2s05_D;
    GUESS(2).I2=I2s05_D;
    GUESS(2).R2=R2s05_D;
    GUESS(2).N2=N2s05_D;
    GUESS(2).uD1=uD1s05;
    GUESS(2).uD2=uD2s05;
    GUESS(2).D1=D1s05;
    GUESS(2).D2=D2s05;
end

%SCAR=0.075;

         
%[Results, solution, ts10_D, S1s10_D, S2s10_D, E1s10_D, E2s10_D, I1s10_D, I2s10_D, R1s10_D, R2s10_D, N1s10_D, N2s10_D, uD1s10, uD2s10, D1s10, D2s10] = ...
 %   covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


if false %Uses results from SCAR=0.1 for a guess for SCAR=0.05
    GUESS(2).S1=S1s10_D;
    GUESS(2).E1=E1s10_D;
    GUESS(2).I1=I1s10_D;
    GUESS(2).R1=R1s10_D;
    GUESS(2).N1=N1s10_D;
    GUESS(2).S2=S2s10_D;
    GUESS(2).E2=E2s10_D;
    GUESS(2).I2=I2s10_D;
    GUESS(2).R2=R2s10_D;
    GUESS(2).N2=N2s10_D;
    GUESS(2).uD1=uD1s10;
    GUESS(2).uD2=uD2s10;
    GUESS(2).D1=D1s10;
    GUESS(2).D2=D2s10;
end

SCAR=0.10;

         
[Results, solution, ts10_D, S1s10_D, S2s10_D, E1s10_D, E2s10_D, I1s10_D, I2s10_D, R1s10_D, R2s10_D, N1s10_D, N2s10_D, uD1s10, uD2s10, D1s10, D2s10] = ...
    covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


if true %Uses results from SCAR=0.1 for a guess for SCAR=0.05
    GUESS(2).S1=S1s10_D;
    GUESS(2).E1=E1s10_D;
    GUESS(2).I1=I1s10_D;
    GUESS(2).R1=R1s10_D;
    GUESS(2).N1=N1s10_D;
    GUESS(2).S2=S2s10_D;
    GUESS(2).E2=E2s10_D;
    GUESS(2).I2=I2s10_D;
    GUESS(2).R2=R2s10_D;
    GUESS(2).N2=N2s10_D;
    GUESS(2).uD1=uD1s10;
    GUESS(2).uD2=uD2s10;
    GUESS(2).D1=D1s10;
    GUESS(2).D2=D2s10;
end




 SCAR=0.14;
    
[Results, solution, ts15_D, S1s15_D, S2s15_D, E1s15_D, E2s15_D, I1s15_D, I2s15_D, R1s15_D, R2s15_D, N1s15_D, N2s15_D, uD1s15, uD2s15, D1s15, D2s15] = ...
    covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS); 


if true %Uses results from SCAR=0.1 for a guess for SCAR=0.15
    GUESS(2).S1=S1s15_D;
    GUESS(2).E1=E1s15_D;
    GUESS(2).I1=I1s15_D;
    GUESS(2).R1=R1s15_D;
    GUESS(2).N1=N1s15_D;
    GUESS(2).S2=S2s15_D;
    GUESS(2).E2=E2s15_D;
    GUESS(2).I2=I2s15_D;
    GUESS(2).R2=R2s15_D;
    GUESS(2).N2=N2s15_D;
    GUESS(2).uD1=uD1s15;
    GUESS(2).uD2=uD2s15;
    GUESS(2).D1=D1s15;
    GUESS(2).D2=D2s15;
end

 SCAR=0.15;
    
[Results, solution, ts15_D, S1s15_D, S2s15_D, E1s15_D, E2s15_D, I1s15_D, I2s15_D, R1s15_D, R2s15_D, N1s15_D, N2s15_D, uD1s15, uD2s15, D1s15, D2s15] = ...
    covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS); 


if false %Uses results from SCAR=0.1 for a guess for SCAR=0.15
    GUESS(2).S1=S1s15_D;
    GUESS(2).E1=E1s15_D;
    GUESS(2).I1=I1s15_D;
    GUESS(2).R1=R1s15_D;
    GUESS(2).N1=N1s15_D;
    GUESS(2).S2=S2s15_D;
    GUESS(2).E2=E2s15_D;
    GUESS(2).I2=I2s15_D;
    GUESS(2).R2=R2s15_D;
    GUESS(2).N2=N2s15_D;
    GUESS(2).uD1=uD1s15;
    GUESS(2).uD2=uD2s15;
    GUESS(2).D1=D1s15;
    GUESS(2).D2=D2s15;
end







end

%%% Allocation of vaccines under different vaccine scarcity
if true
    
    %Initial guess to find unconstained amount
GUESS=[]; % = empty : Initial guess is initial conditions
 
 CASE=3; %=3 Optimal Vaccine
 
   SCAR=0.025;

[Results, solution, ts05, S1s05, S2s05, E1s05, E2s05, I1s05, I2s05, R1s05, R2s05, N1s05, N2s05, uV1s05, uV2s05, V1s05, V2s05] = ...
    covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


  if true %Uses results from SCAR=0.1 for a guess for SCAR=0.05
    GUESS(3).S1=S1s05;
    GUESS(3).E1=E1s05;
    GUESS(3).I1=I1s05;
    GUESS(3).R1=R1s05;
    GUESS(3).N1=N1s05;
    GUESS(3).S2=S2s05;
    GUESS(3).E2=E2s05;
    GUESS(3).I2=I2s05;
    GUESS(3).R2=R2s05;
    GUESS(3).N2=N2s05;
    GUESS(3).uV1=uV1s05;
    GUESS(3).uV2=uV2s05;
    GUESS(3).V1=V1s05;
    GUESS(3).V2=V2s05;
  end
  
 
 
  SCAR=0.05;

[Results, solution, ts05, S1s05, S2s05, E1s05, E2s05, I1s05, I2s05, R1s05, R2s05, N1s05, N2s05, uV1s05, uV2s05, V1s05, V2s05] = ...
    covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


  if true %Uses results from SCAR=0.1 for a guess for SCAR=0.05
    GUESS(3).S1=S1s05;
    GUESS(3).E1=E1s05;
    GUESS(3).I1=I1s05;
    GUESS(3).R1=R1s05;
    GUESS(3).N1=N1s05;
    GUESS(3).S2=S2s05;
    GUESS(3).E2=E2s05;
    GUESS(3).I2=I2s05;
    GUESS(3).R2=R2s05;
    GUESS(3).N2=N2s05;
    GUESS(3).uV1=uV1s05;
    GUESS(3).uV2=uV2s05;
    GUESS(3).V1=V1s05;
    GUESS(3).V2=V2s05;
  end

  

  SCAR=0.1;

[Results, solution, ts10, S1s10, S2s10, E1s10, E2s10, I1s10, I2s10, R1s10, R2s10, N1s10, N2s10, uV1s10, uV2s10, V1s10, V2s10] = ...
    covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

if true %Uses results from SCAR=0.1 for a guess for SCAR=0.05
    GUESS(3).S1=S1s10;
    GUESS(3).E1=E1s10;
    GUESS(3).I1=I1s10;
    GUESS(3).R1=R1s10;
    GUESS(3).N1=N1s10;
    GUESS(3).S2=S2s10;
    GUESS(3).E2=E2s10;
    GUESS(3).I2=I2s10;
    GUESS(3).R2=R2s10;
    GUESS(3).N2=N2s10;
    GUESS(3).uV1=uV1s10;
    GUESS(3).uV2=uV2s10;
    GUESS(3).V1=V1s10;
    GUESS(3).V2=V2s10 ;
end


 
  SCAR=0.15;
 
[Results, solution, ts15, S1s15, S2s15, E1s15, E2s15, I1s15, I2s15, R1s15, R2s15, N1s15, N2s15, uV1s15, uV2s15, V1s15, V2s15] = ...
    covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


if false %Uses results from SCAR=0.1 for a guess for SCAR=0.05
    GUESS(3).S1=S1s15;
    GUESS(3).E1=E1s15;
    GUESS(3).I1=I1s15;
    GUESS(3).R1=R1s15;
    GUESS(3).N1=N1s15;
    GUESS(3).S2=S2s15;
    GUESS(3).E2=E2s15;
    GUESS(3).I2=I2s15;
    GUESS(3).R2=R2s15;
    GUESS(3).N2=N2s15;
    GUESS(3).uV1=uV1s15;
    GUESS(3).uV2=uV2s15;
    GUESS(3).V1=V1s15;
    GUESS(3).V2=V2s15;
end








end

%%% Ad Hoc Allocation of drugs under different drug scarcity
if true
    
%Initial guess to find unconstained amount
GUESS=[]; % = empty : Initial guess is initial conditions

CASE=4; %4= Ad Hoc Symptomatic Drug

SCAR=0.025;
 
    
[Results, solution, ts05_Dah, S1s05_Dah, S2s05_Dah, E1s05_Dah, E2s05_Dah, I1s05_Dah, I2s05_Dah, R1s05_Dah, R2s05_Dah, N1s05_Dah, N2s05_Dah, uD1s05_ah, uD2s05_ah, D1s05_ah, D2s05_ah] = ...
    covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


if true %Uses results from SCAR=0.1 for a guess for SCAR=0.15
    GUESS(2).S1=S1s05_Dah;
    GUESS(2).E1=E1s05_Dah;
    GUESS(2).I1=I1s05_Dah;
    GUESS(2).R1=R1s05_Dah;
    GUESS(2).N1=N1s05_Dah;
    GUESS(2).S2=S2s05_Dah;
    GUESS(2).E2=E2s05_Dah;
    GUESS(2).I2=I2s05_Dah;
    GUESS(2).R2=R2s05_Dah;
    GUESS(2).N2=N2s05_Dah;
    GUESS(2).uD1=uD1s05_ah*2;
    GUESS(2).uD2=uD2s05_ah*2;
    GUESS(2).D1=D1s05_ah;
    GUESS(2).D2=D2s05_ah;
    
    GUESS(4)=GUESS(2);
end


 SCAR=0.05;
 
[Results, solution, ts05_Dah, S1s05_Dah, S2s05_Dah, E1s05_Dah, E2s05_Dah, I1s05_Dah, I2s05_Dah, R1s05_Dah, R2s05_Dah, N1s05_Dah, N2s05_Dah, uD1s05_ah, uD2s05_ah, D1s05_ah, D2s05_ah] = ...
    covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


if true %Uses results from SCAR=0.1 for a guess for SCAR=0.15
    GUESS(2).S1=S1s05_Dah;
    GUESS(2).E1=E1s05_Dah;
    GUESS(2).I1=I1s05_Dah;
    GUESS(2).R1=R1s05_Dah;
    GUESS(2).N1=N1s05_Dah;
    GUESS(2).S2=S2s05_Dah;
    GUESS(2).E2=E2s05_Dah;
    GUESS(2).I2=I2s05_Dah;
    GUESS(2).R2=R2s05_Dah;
    GUESS(2).N2=N2s05_Dah;
    GUESS(2).uD1=uD1s05_ah*2;
    GUESS(2).uD2=uD2s05_ah*2;
    GUESS(2).D1=D1s05_ah;
    GUESS(2).D2=D2s05_ah;
    
    GUESS(4)=GUESS(2);
end


 SCAR=0.1;
    
[Results, solution, ts10_Dah, S1s10_Dah, S2s10_Dah, E1s10_Dah, E2s10_Dah, I1s10_Dah, I2s10_Dah, R1s10_Dah, R2s10_Dah, N1s10_Dah, N2s10_Dah, uD1s10_ah, uD2s10_ah, D1s10_ah, D2s10_ah] = ...
    covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


if true %Uses results from SCAR=0.1 for a guess for SCAR=0.05
    GUESS(2).S1=S1s10_Dah;
    GUESS(2).E1=E1s10_Dah;
    GUESS(2).I1=I1s10_Dah;
    GUESS(2).R1=R1s10_Dah;
    GUESS(2).N1=N1s10_Dah;
    GUESS(2).S2=S2s10_Dah;
    GUESS(2).E2=E2s10_Dah;
    GUESS(2).I2=I2s10_Dah;
    GUESS(2).R2=R2s10_Dah;
    GUESS(2).N2=N2s10_Dah;
    GUESS(2).uD1=uD1s10_ah;
    GUESS(2).uD2=uD2s10_ah;
    GUESS(2).D1=D1s10_ah;
    GUESS(2).D2=D2s10_ah;
    
    GUESS(4)=GUESS(2);
end



 SCAR=0.15;
    
[Results, solution, ts15_Dah, S1s15_Dah, S2s15_Dah, E1s15_Dah, E2s15_Dah, I1s15_Dah, I2s15_Dah, R1s15_Dah, R2s15_Dah, N1s15_Dah, N2s15_Dah, uD1s15_ah, uD2s15_ah, D1s15_ah, D2s01_ah] = ...
    covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);
 

end

%%% Ad Hoc Allocation of vaccines under different vaccine scarcity
if true
    
 %Initial guess to find unconstained amount
GUESS=[]; % = empty : Initial guess is initial conditions

 CASE=5; %=5 Ad Hoc Vaccine
 
   SCAR=0.025;
    
[Results_ah, solution_ah, ts05_ah, S1s05_ah, S2s05_ah, E1s05_ah, E2s05_ah, I1s05_ah, I2s05_ah, R1s05_ah, R2s05_ah, N1s05_ah, N2s05_ah, uV1s05_ah, uV2s05_ah, V1s05_ah, V2s05_ah] = ...
    covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

   if true %Uses results from SCAR=0.1 for a guess for SCAR=0.05
    GUESS(3).S1=S1s05_ah;
    GUESS(3).E1=E1s05_ah;
    GUESS(3).I1=I1s05_ah;
    GUESS(3).R1=R1s05_ah;
    GUESS(3).N1=N1s05_ah;
    GUESS(3).S2=S2s05_ah;
    GUESS(3).E2=E2s05_ah;
    GUESS(3).I2=I2s05_ah;
    GUESS(3).R2=R2s05_ah;
    GUESS(3).N2=N2s05_ah;
    GUESS(3).uV1=uV1s05_ah;
    GUESS(3).uV2=uV2s05_ah;
    GUESS(3).V1=V1s05_ah;
    GUESS(3).V2=V2s05_ah;
    
    GUESS(5)=GUESS(3);
  end
  
 
  SCAR=0.05;
    
[Results_ah, solution_ah, ts05_ah, S1s05_ah, S2s05_ah, E1s05_ah, E2s05_ah, I1s05_ah, I2s05_ah, R1s05_ah, R2s05_ah, N1s05_ah, N2s05_ah, uV1s05_ah, uV2s05_ah, V1s05_ah, V2s05_ah] = ...
    covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

 
   if true %Uses results from SCAR=0.1 for a guess for SCAR=0.05
    GUESS(3).S1=S1s05_ah;
    GUESS(3).E1=E1s05_ah;
    GUESS(3).I1=I1s05_ah;
    GUESS(3).R1=R1s05_ah;
    GUESS(3).N1=N1s05_ah;
    GUESS(3).S2=S2s05_ah;
    GUESS(3).E2=E2s05_ah;
    GUESS(3).I2=I2s05_ah;
    GUESS(3).R2=R2s05_ah;
    GUESS(3).N2=N2s05_ah;
    GUESS(3).uV1=uV1s05_ah;
    GUESS(3).uV2=uV2s05_ah;
    GUESS(3).V1=V1s05_ah;
    GUESS(3).V2=V2s05_ah;
    
    GUESS(5)=GUESS(3);
  
   end
  

 SCAR=0.1;

[Results_ah, solution_ah, ts10_ah, S1s10_ah, S2s10_ah, E1s10_ah, E2s10_ah, I1s10_ah, I2s10_ah, R1s10_ah, R2s10_ah, N1s10_ah, N2s10_ah, uV1s10_ah, uV2s10_ah, V1s10_ah, V2s10_ah] = ...
    covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


if true %Uses results from SCAR=0.1 for a guess for SCAR=0.05
    GUESS(3).S1=S1s10_ah;
    GUESS(3).E1=E1s10_ah;
    GUESS(3).I1=I1s10_ah;
    GUESS(3).R1=R1s10_ah;
    GUESS(3).N1=N1s10_ah;
    GUESS(3).S2=S2s10_ah;
    GUESS(3).E2=E2s10_ah;
    GUESS(3).I2=I2s10_ah;
    GUESS(3).R2=R2s10_ah;
    GUESS(3).N2=N2s10_ah;
    GUESS(3).uV1=uV1s10_ah;
    GUESS(3).uV2=uV2s10_ah;
    GUESS(3).V1=V1s10_ah;
    GUESS(3).V2=V2s10_ah ;
    
    GUESS(5)=GUESS(3);

end



  SCAR=0.15;

[Results_ah, solution_ah, ts15_ah, S1s15_ah, S2s15_ah, E1s15_ah, E2s15_ah, I1s15_ah, I2s15_ah, R1s15_ah, R2s15_ah, N1s15_ah, N2s15_ah, uV1s15_ah, uV2s15_ah, V1s15_ah, V2s01_ah] = ...
    covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%   Sensitivity Analyses   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SA Adjustment cost parameter with Vaccine
if false
    
    [N1,N2,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,r,cD,cV,cAdj,cI,w,params,R0] = covid_Parameters_months()


    ODE=1;
    omega=0;%1/6;
    SCAR=0.1;
    
    CADJ=[10e00 10e01 10e02 10e03 10e04 10e05 10e06 10e07 10e08];
    

    %Initial guess to find unconstained amount
    GUESS=[]; % = empty : Initial guess is initial conditions
    
    
VaccineProp=zeros(1,length(CADJ));
VaccineInfection=zeros(1,length(CADJ));
CostProp=zeros(1,length(CADJ));

    for i=1:length(CADJ)
        
        cAdj=CADJ(i);

CASE=3;
    
[Results, solution, ts10, S1s10, S2s10, E1s10, E2s10, I1s10, I2s10, R1s10, R2s10, N1s10, N2s10, uV1s10, uV2s10, V1s10, V2s10] = ...
    covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


CASE=5;

[Results_ah, solution_ah, ts10_ah, S1s10_ah, S2s10_ah, E1s10_ah, E2s10_ah, I1s10_ah, I2s10_ah, R1s10_ah, R2s10_ah, N1s10_ah, N2s10_ah, uV1s10_ah, uV2s10_ah, V1s10_ah, V2s10_ah] = ...
    covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


%Type of interpolation
fitApprox = fittype('pchipinterp'); %pchipinterp

ts=ts10;

%Interpolation of Cumulative cases
gI_Optimal=fit(ts,(I1s10+I2s10),fitApprox);
gI_AdHoc=fit(ts,(I1s10_ah+I2s10_ah),fitApprox);
CumI_Optimal=integrate(gI_Optimal,T,0);
CumI_AdHoc=integrate(gI_AdHoc,T,0);

%Terminal time period where treatment becomes nonbinding
scarce=max(find(round(uV1s10+uV2s10,4)==round(MaxTreat(2)*SCAR,4)));

%VaccineProp(i)=CumV1_Optimal./CumV2_Optimal
VaccineProp(i)=var((uV1s10(1:scarce)-uV1s10_ah(1:scarce))./uV1s10_ah(1:scarce))
VaccineInfection(i)=(CumI_Optimal-CumI_AdHoc)./CumI_AdHoc

%Interpolation of Cost
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cAdj.*((N2s10./(N1s10+N2s10)).*uV1s10 - (N1s10./(N1s10+N2s10)).*uV2s10).^2)),fitApprox); 
gVaccineCost_Optimal = fit(ts,(exp(-r.*ts).*(cV.*(uV1s10+uV2s10))),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);
VaccineCost_Optimal=integrate(gVaccineCost_Optimal,T,0);

CostProp(i)=Penalty_Optimal/VaccineCost_Optimal

    end
    
    
    fig=figure
    subplot(211)
    yyaxis left
    plot(log(CADJ),VaccineProp*100,'LineWidth',3);hold on
    xticks(log(CADJ));
    xticklabels({'10','100','1K','Base Case','100K','1M','VSL','100M', '1B'});
    xlim([log(CADJ(1)) log(CADJ(end))])
    ylabel({'Variance of the Optimal ','Deviation From the Ad Hoc'},'Interpreter','tex','FontSize', 12);     
    plot([log(CADJ(4)) log(CADJ(4))],[VaccineProp(4)*100 VaccineProp(4)*100],'k.','LineWidth',1,'MarkerSize',24); 
    
    yyaxis right
    plot(log(CADJ),VaccineInfection*100,'LineWidth',3);hold on
    ylabel({'Difference in','Cumulative Cases'},'Interpreter','tex','FontSize', 12);
    plot([log(CADJ(4)) log(CADJ(4))],[VaccineInfection(4)*100 VaccineInfection(4)*100],'k.','LineWidth',1,'MarkerSize',24); 
    
    subplot(212)
    plot(log(CADJ),CostProp,'LineWidth',3);hold on
    xticks(log(CADJ));
    xticklabels({'10','100','1K','Base Case','100K','1M','VSL','100M', '1B'});
    xlim([log(CADJ(1)) log(CADJ(end))])
    plot([log(CADJ(4)) log(CADJ(4))],[CostProp(4) CostProp(4)],'k.','LineWidth',1,'MarkerSize',24);
    ylabel({'$\frac{\textrm{Total Workability}}{\textrm{Total Vaccine}}$'},'Interpreter','latex', 'FontSize', 12);
    
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Workability Cost Parameter', 'FontSize', 16,'Interpreter','tex');
    
    
    %saveas(gcf,'SA_WorkabilityCost.png'); hold off
end

% SA Adjustment cost parameter with Drugs
if false
    
    [N1,N2,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,r,cD,cV,cAdj,cI,w,params,R0] = covid_Parameters_months()


    ODE=2;%2;
    omega=1/6;
    SCAR=0.1;
    
    CADJ=[10e00 10e01 10e02 10e03 10e05 10e06 10e07 10e08];
    

    %Initial guess to find unconstained amount
    GUESS=[]; % = empty : Initial guess is initial conditions
    
    
DrugProp=zeros(1,length(CADJ));
DrugProp2=zeros(1,length(CADJ));
DrugInfection=zeros(1,length(CADJ));
CostProp=zeros(1,length(CADJ));

    for i=1:length(CADJ)
        
        cAdj=CADJ(i);
 
if ODE==2
       if i==2 || i==3 || i==4 || i==5 || i==6 || i==7 || i==8
           
           if true 
    GUESS(2).S1=S1s10_D;
    GUESS(2).E1=E1s10_D;
    GUESS(2).I1=I1s10_D;
    GUESS(2).R1=R1s10_D;
    GUESS(2).N1=N1s10_D;
    GUESS(2).S2=S2s10_D;
    GUESS(2).E2=E2s10_D;
    GUESS(2).I2=I2s10_D;
    GUESS(2).R2=R2s10_D;
    GUESS(2).N2=N2s10_D;
    GUESS(2).uD1=uD1s10;
    GUESS(2).uD2=uD2s10;
    GUESS(2).D1=D1s10;
    GUESS(2).D2=D2s10;
    
    GUESS(4)=GUESS(2);
           end

       end
end
CASE=2;

if ODE==2
SCAR=0.09;

if omega==1/6
SCAR=0.1;
end

[Results, solution, ts10_D, S1s10_D, S2s10_D, E1s10_D, E2s10_D, I1s10_D, I2s10_D, R1s10_D, R2s10_D, N1s10_D, N2s10_D, uD1s10, uD2s10, D1s10, D2s10] = ...
    covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

if true %Uses results from SCAR=0.09 for a guess for SCAR=0.1
    GUESS(2).S1=S1s10_D;
    GUESS(2).E1=E1s10_D;
    GUESS(2).I1=I1s10_D;
    GUESS(2).R1=R1s10_D;
    GUESS(2).N1=N1s10_D;
    GUESS(2).S2=S2s10_D;
    GUESS(2).E2=E2s10_D;
    GUESS(2).I2=I2s10_D;
    GUESS(2).R2=R2s10_D;
    GUESS(2).N2=N2s10_D;
    GUESS(2).uD1=uD1s10;
    GUESS(2).uD2=uD2s10;
    GUESS(2).D1=D1s10;
    GUESS(2).D2=D2s10;
    
    GUESS(4)=GUESS(2);
end
end 


[Results, solution, ts10_D, S1s10_D, S2s10_D, E1s10_D, E2s10_D, I1s10_D, I2s10_D, R1s10_D, R2s10_D, N1s10_D, N2s10_D, uD1s10, uD2s10, D1s10, D2s10] = ...
    covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


CASE=4;

[Results, solution, ts10_Dah, S1s10_Dah, S2s10_Dah, E1s10_Dah, E2s10_Dah, I1s10_Dah, I2s10_Dah, R1s10_Dah, R2s10_Dah, N1s10_Dah, N2s10_Dah, uD1s10_ah, uD2s10_ah, D1s10_ah, D2s10_ah] = ...
    covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


%Type of interpolation
fitApprox = fittype('pchipinterp'); %pchipinterp

ts=ts10_D;

%Interpolation of Cumulative cases
gI_Optimal=fit(ts,(I1s10_D+I2s10_D),fitApprox);
gI_AdHoc=fit(ts,(I1s10_Dah+I2s10_Dah),fitApprox);
CumI_Optimal=integrate(gI_Optimal,T,0);
CumI_AdHoc=integrate(gI_AdHoc,T,0);

%Terminal time period where treatment becomes nonbinding
scarce=max(find(round(uD1s10+uD2s10,4)==round(MaxTreat(1)*SCAR,4)));

%VaccineProp(i)=CumV1_Optimal./CumV2_Optimal
DrugProp(i)=var((uD1s10(1:scarce)-uD1s10_ah(1:scarce)))%./uD1s10_ah(1:scarce))
DrugProp2(i)=var((uD2s10(1:scarce)-uD2s10_ah(1:scarce)))%./uD2s10_ah(1:scarce))
DrugInfection(i)=(CumI_Optimal-CumI_AdHoc)./CumI_AdHoc

%Interpolation of Cost
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cAdj.*((I2s10_D./(I1s10_D+I2s10_D)).*uD1s10 - (I1s10_D./(I1s10_D+I2s10_D)).*uD2s10).^2)),fitApprox); 
gDrugCost_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(uD1s10+uD2s10))),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);
DrugCost_Optimal=integrate(gDrugCost_Optimal,T,0);

CostProp(i)=Penalty_Optimal/DrugCost_Optimal

    end
    
    
    fig=figure
    subplot(211)
    yyaxis left
    plot(log(CADJ),DrugProp*100,'LineWidth',3);hold on
    xticks(log(CADJ));
    xticklabels({'10','100','1K','Base Case','100K','1M','VSL','100M', '1B'});
    xlim([log(CADJ(1)) log(CADJ(end))])
    ylabel({'Variance of the Optimal ','Deviation From the Ad Hoc'},'Interpreter','tex','FontSize', 12);     
    plot([log(CADJ(4)) log(CADJ(4))],[DrugProp(4)*100 DrugProp(4)*100],'k.','LineWidth',1,'MarkerSize',24); 
    
    yyaxis right
    plot(log(CADJ),DrugInfection*100,'LineWidth',3);hold on
    ylabel({'Difference in','Cumulative Cases'},'Interpreter','tex','FontSize', 12);
    plot([log(CADJ(4)) log(CADJ(4))],[DrugInfection(4)*100 DrugInfection(4)*100],'k.','LineWidth',1,'MarkerSize',24); 
    
    subplot(212)
    plot(log(CADJ),CostProp,'LineWidth',3);hold on
    xticks(log(CADJ));
    xticklabels({'10','100','1K','Base Case','100K','1M','VSL','100M', '1B'});
    xlim([log(CADJ(1)) log(CADJ(end))])
    plot([log(CADJ(4)) log(CADJ(4))],[CostProp(4) CostProp(4)],'k.','LineWidth',1,'MarkerSize',24);
    ylabel({'$\frac{\textrm{Total Workability}}{\textrm{Total Vaccine}}$'},'Interpreter','latex', 'FontSize', 12);
    
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Workability Cost Parameter', 'FontSize', 16,'Interpreter','tex');
    
    
    %saveas(gcf,'SA_WorkabilityCost.png'); hold off
end



%Vaccine effectiveness
if false
    
        %[N1,N2,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,r,cD,cV,cAdj,cI,w,params,R0] = covid_Parameters_months()

        
   
    ODE=2;
    omega=1/6;
    SCAR=0.1;

    
    QV=[0.55 0.65 .75 .85 .95];
    

    %Initial guess to find unconstained amount
    GUESS=[]; % = empty : Initial guess is initial conditions
    
    
VaccineProp=zeros(1,length(QV));
VaccineInfection=zeros(1,length(QV));
CostProp=zeros(1,length(QV));

    for i=1:length(QV)
        
        qV=QV(i);
        
CASE=3;
    
[Results, solution, ts10, S1s10, S2s10, E1s10, E2s10, I1s10, I2s10, R1s10, R2s10, N1s10, N2s10, uV1s10, uV2s10, V1s10, V2s10] = ...
    covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


CASE=5;

[Results_ah, solution_ah, ts10_ah, S1s10_ah, S2s10_ah, E1s10_ah, E2s10_ah, I1s10_ah, I2s10_ah, R1s10_ah, R2s10_ah, N1s10_ah, N2s10_ah, uV1s10_ah, uV2s10_ah, V1s10_ah, V2s10_ah] = ...
    covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


%Type of interpolation
fitApprox = fittype('pchipinterp'); %pchipinterp

ts=ts10;
%Interpolation of Cumulative cases
gI_Optimal=fit(ts,(I1s10+I2s10),fitApprox);
gI_AdHoc=fit(ts,(I1s10_ah+I2s10_ah),fitApprox);
CumI_Optimal=integrate(gI_Optimal,T,0);
CumI_AdHoc=integrate(gI_AdHoc,T,0);

%Terminal time period where treatment becomes nonbinding
scarce=max(find(round(uV1s10+uV2s10,4)==round(MaxTreat(2)*SCAR,4)));

%VaccineProp(i)=CumV1_Optimal./CumV2_Optimal
VaccineProp(i)=var((uV1s10(1:scarce)-uV1s10_ah(1:scarce))./uV1s10_ah(1:scarce))
VaccineInfection(i)=(CumI_Optimal-CumI_AdHoc)./CumI_AdHoc

%Interpolation of Cost
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cAdj.*((N2s10./(N1s10+N2s10)).*uV1s10 - (N1s10./(N1s10+N2s10)).*uV2s10).^2)),fitApprox); 
gVaccineCost_Optimal = fit(ts,(exp(-r.*ts).*(cV.*(uV1s10+uV2s10))),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);
VaccineCost_Optimal=integrate(gVaccineCost_Optimal,T,0);

CostProp(i)=Penalty_Optimal/VaccineCost_Optimal

    end
    
    
   fig=figure
    subplot(211)
    yyaxis left
    plot(QV,VaccineProp,'LineWidth',3);hold on
    xticks(QV);
    xticklabels({'0.55','Base Case', '0.75', '0.85', '0.95'});
    xlim([QV(1) QV(end)])
%    yticks([1 1.1 VaccineProp(2) 1.4]);
 %   yticklabels({'1','1.1','1.22','1.4'});
      
 %   plot([QV(2) QV(2)],[VaccineProp(2) VaccineProp(2)],'k.','LineWidth',1,'MarkerSize',35);
  %  plot([0 QV(2)],[VaccineProp(2) VaccineProp(2)],'k','LineWidth',1);
  %  plot([QV(2) QV(2)],[1 VaccineProp(2)],'k','LineWidth',1);
    
    yyaxis right
    plot(QV,VaccineInfection,'LineWidth',3);hold on
    
    title({'Proportion of Total Vaccines in State 1','Total Vaccines in State 2'}, 'FontSize', 16);
      
    subplot(212)
    plot(QV,CostProp,'LineWidth',3);hold on
    xticks(QV);
    xticklabels({'0.55','Base Case', '0.75', '0.85', '0.95'});
    xlim([QV(1) QV(end)])
    yticks([0 CostProp(2) 3 5]);
    yticklabels({'0','1.3','3','5'});
    
    plot([QV(2) QV(2)],[CostProp(2) CostProp(2)],'k.','LineWidth',1,'MarkerSize',35);
    plot([0 QV(2)],[CostProp(2) CostProp(2)],'k','LineWidth',1); hold on
    plot([QV(2) QV(2)],[0 CostProp(2)],'k','LineWidth',1);
    title({'Proportion of Total Workability Costs','to Total Vaccine Costs'}, 'FontSize', 16);
    
     han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Effectiveness of Vaccine', 'FontSize', 16);
       saveas(gcf,'SA_EffectivenessVaccine.png'); hold off
end

%Drug effectiveness
if false
    
        
   
    ODE=1;
    omega=1/6;
    SCAR=0.1;
    %Note that these results were obtained using "fem1s" collocation points
    %instead of Gauss; this need to be changed manually in the "covid19_Drugs.m" file
    
    if ODE==1 && omega==0
    Nset=90;
    end
   QD=[.55 0.65 0.75 .85 .95];


    %Initial guess to find unconstained amount
    GUESS=[]; % = empty : Initial guess is initial conditions
    
    
DrugProp=zeros(1,length(QD));
DrugProp2=zeros(1,length(QD));
DrugInfection=zeros(1,length(QD));
CostProp=zeros(1,length(QD));

Penalty_Opt=zeros(1,length(QD));
DrugCost_Opt=zeros(1,length(QD));


DrugPaths=zeros(Nset+2,length(QD));
DrugPaths2=zeros(Nset+2,length(QD));

DrugPathsAH=zeros(Nset+2,length(QD));
DrugPaths2AH=zeros(Nset+2,length(QD));

    for i=1:length(QD)
        
        qD=QD(i);
        
if ODE==1
       if i==2 || i==3 || i==4 || i==5 || i==6 || i==7 || i==8 || i==9
           
           if true 
    GUESS(2).S1=S1s10_D;
    GUESS(2).E1=E1s10_D;
    GUESS(2).I1=I1s10_D;
    GUESS(2).R1=R1s10_D;
    GUESS(2).N1=N1s10_D;
    GUESS(2).S2=S2s10_D;
    GUESS(2).E2=E2s10_D;
    GUESS(2).I2=I2s10_D;
    GUESS(2).R2=R2s10_D;
    GUESS(2).N2=N2s10_D;
    GUESS(2).uD1=uD1s10;
    GUESS(2).uD2=uD2s10;
    GUESS(2).D1=D1s10;
    GUESS(2).D2=D2s10;
    
    GUESS(4)=GUESS(2);
           end

       end
end
CASE=2;

if ODE==1
SCAR=0.09;
 
[Results, solution, ts10_D, S1s10_D, S2s10_D, E1s10_D, E2s10_D, I1s10_D, I2s10_D, R1s10_D, R2s10_D, N1s10_D, N2s10_D, uD1s10, uD2s10, D1s10, D2s10] = ...
    covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

if true %Uses results from SCAR=0.09 for a guess for SCAR=0.1
    GUESS(2).S1=S1s10_D;
    GUESS(2).E1=E1s10_D;
    GUESS(2).I1=I1s10_D;
    GUESS(2).R1=R1s10_D;
    GUESS(2).N1=N1s10_D;
    GUESS(2).S2=S2s10_D;
    GUESS(2).E2=E2s10_D;
    GUESS(2).I2=I2s10_D;
    GUESS(2).R2=R2s10_D;
    GUESS(2).N2=N2s10_D;
    GUESS(2).uD1=uD1s10;
    GUESS(2).uD2=uD2s10;
    GUESS(2).D1=D1s10;
    GUESS(2).D2=D2s10;
    
    GUESS(4)=GUESS(2);
end
end 

SCAR=0.1;

         
[Results, solution, ts10_D, S1s10_D, S2s10_D, E1s10_D, E2s10_D, I1s10_D, I2s10_D, R1s10_D, R2s10_D, N1s10_D, N2s10_D, uD1s10, uD2s10, D1s10, D2s10] = ...
    covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


CASE=4;
    
[Results, solution, ts10_Dah, S1s10_Dah, S2s10_Dah, E1s10_Dah, E2s10_Dah, I1s10_Dah, I2s10_Dah, R1s10_Dah, R2s10_Dah, N1s10_Dah, N2s10_Dah, uD1s10_ah, uD2s10_ah, D1s10_ah, D2s10_ah] = ...
    covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


%Type of interpolation
fitApprox = fittype('pchipinterp'); %pchipinterp

ts=ts10_D;
%Interpolation of Cumulative cases
gI_Optimal=fit(ts,(I1s10_D+I2s10_D),fitApprox);
gI_AdHoc=fit(ts,(I1s10_Dah+I2s10_Dah),fitApprox);
CumI_Optimal=integrate(gI_Optimal,T,0);
CumI_AdHoc=integrate(gI_AdHoc,T,0);

%Terminal time period where treatment becomes nonbinding
scarce=max(find(round(uD1s10+uD2s10,4)==round(MaxTreat(1)*SCAR,4)));

%VaccineProp(i)=CumV1_Optimal./CumV2_Optimal
DrugProp(i)=var((uD1s10(1:scarce)-uD1s10_ah(1:scarce)))%./uD1s10_ah(1:scarce))
DrugProp2(i)=var((uD2s10(1:scarce)-uD2s10_ah(1:scarce)))%./uD2s10_ah(1:scarce))
DrugInfection(i)=(CumI_Optimal-CumI_AdHoc)./CumI_AdHoc

%
DrugPaths(:,i)=uD1s10;
DrugPaths2(:,i)=uD2s10;
DrugPathsAH(:,i)=uD1s10_ah;
DrugPaths2AH(:,i)=uD2s10_ah;


%Interpolation of Cost
gPenalty_Optimal = fit(ts,(exp(-r.*ts).*(cAdj.*((I2s10_D./(I1s10_D+I2s10_D)).*uD1s10 - (I1s10_D./(I1s10_D+I2s10_D)).*uD2s10).^2)),fitApprox); 
gDrugCost_Optimal = fit(ts,(exp(-r.*ts).*(cD.*(uD1s10+uD2s10))),fitApprox); 
Penalty_Optimal=integrate(gPenalty_Optimal,T,0);
DrugCost_Optimal=integrate(gDrugCost_Optimal,T,0);

Penalty_Opt(i)=Penalty_Optimal;
DrugCost_Opt(i)=DrugCost_Optimal;

CostProp(i)=Penalty_Optimal/DrugCost_Optimal

    end
    
    
   fig=figure
    subplot(211)
    yyaxis left
    plot(QD,DrugProp,'LineWidth',3);hold on
    xticks(QD);
    xticklabels({'0.55','Base Case', '0.75', '0.85', '0.95'});
    xlim([QD(1) QD(end)])
%    yticks([1 1.1 VaccineProp(2) 1.4]);
 %   yticklabels({'1','1.1','1.22','1.4'});
      
 %   plot([QV(2) QV(2)],[VaccineProp(2) VaccineProp(2)],'k.','LineWidth',1,'MarkerSize',35);
  %  plot([0 QV(2)],[VaccineProp(2) VaccineProp(2)],'k','LineWidth',1);
  %  plot([QV(2) QV(2)],[1 VaccineProp(2)],'k','LineWidth',1);
    
    yyaxis right
    plot(QD,DrugInfection,'LineWidth',3);hold on
    
    title({'Proportion of Total Vaccines in State 1','Total Vaccines in State 2'}, 'FontSize', 16);
      
    subplot(212)
    plot(QD,CostProp,'LineWidth',3);hold on
    xticks(QD);
    xticklabels({'0.55','Base Case', '0.75', '0.85', '0.95'});
    xlim([QD(1) QD(end)])
    yticks([0 CostProp(2) 3 5]);
    yticklabels({'0','1.3','3','5'});
    
    plot([QD(2) QD(2)],[CostProp(2) CostProp(2)],'k.','LineWidth',1,'MarkerSize',35);
    plot([0 QD(2)],[CostProp(2) CostProp(2)],'k','LineWidth',1); hold on
    plot([QD(2) QD(2)],[0 CostProp(2)],'k','LineWidth',1);
    title({'Proportion of Total Workability Costs','to Total Vaccine Costs'}, 'FontSize', 16);
    
     han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Effectiveness of Vaccine', 'FontSize', 16);
       saveas(gcf,'SA_EffectivenessDrug.png'); hold off
end


