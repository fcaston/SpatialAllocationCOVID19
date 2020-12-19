clear all; close all
%Calling parameter values
[N1,N2,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,r,cD,cV,cAdj,cI,w,params,R0] = covid_Parameters_months()

% Number of time periods
T =4; %(in months)

%Number of collocation points
Nset= 60; 
 
%Objective Function
OBJ=1; %=1 Minimize damages & costs (currently no other obj defined)

%Initial conditions 
%x0ic= [0.9572*N1; 0.0102*N1; 0.0112*N1; 0.0213*N1; 0.9134*N2; 0.0199*N2; 0.0224*N2; 0.0441*N2; N1; N2; 0; 0]; %Simulated using Tian et al (2020)
%x0ic=  [0.9493*N1; 0.0061*N1; 0.0082*N1; 0.0357*N1; 0.9074*N2; 0.0103*N2; 0.0143*N2; 0.0667*N2; N1; N2; 0; 0];
x0ic=  [0.9074*N1; 0.0103*N1; 0.0143*N1; 0.0667*N1; 0.8662*N2; 0.0138*N2; 0.0196*N2; 0.0986*N2; N1; N2; 0; 0];

%Key varying factors
%omega=0  ; ODE=1; %Permanent immunity, compliance to SIP
%omega=0  ; ODE=2; %Permanent immunity, no compliance to SIP
%omega=1/6; ODE=1; %6-month immunity, compliance to SIP
%omega=1/6; ODE=2; %6-month immunity, no compliance to SIP 

%Maximum phyical constraint (does not account for potential scarcity)
%MaxTreat=[0.1.*(N1+N2);N1+N2]; %constraint is technically 0.200908668897642/2 constraint for uncontrolloed R0 of 3.2
%MaxTreat=[0.0317.*(N1+N2);N1+N2]; %constraint for uncontrolled R0 of 2.2
 
MaxTreat=[0.1.*(N1+N2);N1+N2]; %constraint for uncontrolled R0 of 2.2 (it's 0.1091 in fact)


%Evolution of state variables for the no control case
if true
    
CASE=1;
    
SCAR=1;
    
%Initial guess to find unconstained amount
GUESS=[]; % = empty : Initial guess is initial conditions
    
omega=0; ODE=1; %Permanent immunity, compliance to SIP 
    
    [Results, solution, ts, S1s_Dno1PS, S2s_Dno1PS, E1s_Dno1PS, E2s_Dno1PS, I1s_Dno1PS, I2s_Dno1PS, R1s_Dno1PS, R2s_Dno1PS, N1s_Dno1PS, N2s_Dno1PS, uD1s, uD2s, D1s, D2s] = ...
        covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

    [Results, solution, ts, S1s_no1PS, S2s_no1PS, E1s_no1PS, E2s_no1PS, I1s_no1PS, I2s_no1PS, R1s_no1PS, R2s_no1PS, N1s_no1PS, N2s_no1PS, uV1s, uV2s, V1s, V2s] = ...
        covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

omega=0; ODE=2;%6-month immunity, compliance to SIP 
    
    [Results, solution, ts, S1s_Dno3PN, S2s_Dno3PN, E1s_Dno3PN, E2s_Dno3PN, I1s_Dno3PN, I2s_Dno3PN, R1s_Dno3PN, R2s_Dno3PN, N1s_Dno3PN, N2s_Dno3PN, uD1s, uD2s, D1s, D2s] = ...
        covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

    [Results, solution, ts, S1s_no3PN, S2s_no3PN, E1s_no3PN, E2s_no3PN, I1s_no3PN, I2s_no3PN, R1s_no3PN, R2s_no3PN, N1s_no3PN, N2s_no3PN, uV1s, uV2s, V1s, V2s] = ...
        covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);
    
omega=1/6; ODE=1; 
    
    [Results, solution, ts, S1s_Dno2SS, S2s_Dno2SS, E1s_Dno2SS, E2s_Dno2SS, I1s_Dno2SS, I2s_Dno2SS, R1s_Dno2SS, R2s_Dno2SS, N1s_Dno2SS, N2s_Dno2SS, uD1s, uD2s, D1s, D2s] = ...
        covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

    [Results, solution, ts, S1s_no2SS, S2s_no2SS, E1s_no2SS, E2s_no2SS, I1s_no2SS, I2s_no2SS, R1s_no2SS, R2s_no2SS, N1s_no2SS, N2s_no2SS, uV1s, uV2s, V1s, V2s] = ...
        covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);
    
omega=1/6; ODE=2; %6-month immunity, no compliance to SIP
    
        [Results, solution, ts, S1s_Dno4SN, S2s_Dno4SN, E1s_Dno4SN, E2s_Dno4SN, I1s_Dno4SN, I2s_Dno4SN, R1s_Dno4SN, R2s_Dno4SN, N1s_Dno4SN, N2s_Dno4SN, uD1s, uD2s, D1s, D2s] = ...
        covid19_Drugs(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);

    [Results, solution, ts, S1s_no4SN, S2s_no4SN, E1s_no4SN, E2s_no4SN, I1s_no4SN, I2s_no4SN, R1s_no4SN, R2s_no4SN, N1s_no4SN, N2s_no4SN, uV1s, uV2s, V1s, V2s] = ...
        covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,cD,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);
    
    
end


if false
    
    
    
%%% Allocation of symptomatic drugs under different drug scarcity
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
