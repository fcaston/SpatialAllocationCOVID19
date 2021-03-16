function [Results, solution, ts, S1s, S2s, E1s, E2s, I1s, I2s, R1s, R2s, N1s, N2s, uV1s, uV2s, V1s, V2s] = ...
    covid19_Vaccines(T,r,Nset,x0ic,cI,w,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi1,phi2,qV,cV,cAdj,ODE,CASE,MaxTreat,SCAR,OBJ,GUESS);


toms t;
p=tomPhase('p', t, 0, T, Nset,[],'gauss'); %'gauss' 'cheb' or 'fem1s' or 'fem1'
setPhase(p);
tomStates S1 E1 I1 R1 S2 E2 I2 R2 N1 N2 V1 V2
tomControls uV1 uV2

state = [S1; E1; I1; R1; S2; E2; I2; R2; N1; N2; V1; V2];
control = [uV1; uV2];

%Terminal Conditions
cterm=final({});

% Box constraints 
lowerS = zeros(12,1); %Nonnegativity constraints on states
lowerC = zeros(2,1); %Nonnegativity constraints on controls

%ODEs/State equations
if ODE==1 % Compliance to Shelter-in-Place order
    
ode = collocate({
    dot(S1) == omega.*R1 - beta11.*S1.*(I1./N1) - qV.*uV1 ; %S1 dot
    dot(E1) == beta11.*S1.*(I1./N1) - sigma.*E1 ; %E1 dot
    dot(I1) == sigma.*E1 - gamma.*I1 -  phi1.*I1 ; %I1 dot
    dot(R1) == qV.*uV1 + gamma.*I1 - omega.*R1 ; %R1 dot
    dot(S2) == omega.*R2 - beta22.*S2.*(I2./N2) - qV.*uV2 ; %S2 dot
    dot(E2) == beta22.*S2.*(I2./N2) - sigma.*E2 ; %E2 dot
    dot(I2) == sigma.*E2 - gamma.*I2 -  phi2.*I2 ; %I2 dot
    dot(R2) == qV.*uV2  + gamma.*I2 - omega.*R2 ; %R2 dot
    dot(N1) == -phi1.*I1 ;
    dot(N2) == -phi2.*I2 ;
    dot(V1) ==  uV1 ;
    dot(V2) ==  uV2 ;});    
    
end

if ODE==2 % No Compliance to SIP order(Cross contamination)
    
    ode = collocate({
    dot(S1) == omega.*R1 - beta11.*S1.*(I1./N1) - beta12.*S1.*(I2./N2) - qV.*uV1 ; %S1 dot
    dot(E1) == beta11.*S1.*(I1./N1) +  beta12.*S1.*(I2./N2) - sigma.*E1 ; %E1 dot
    dot(I1) == sigma.*E1 - gamma.*I1 -  phi1.*I1 ; %I1 dot
    dot(R1) == qV.*uV1 + gamma.*I1 - omega.*R1 ; %R1 dot
    dot(S2) == omega.*R2 - beta22.*S2.*(I2./N2) - beta21.*S2.*(I1./N1) - qV.*uV2 ; %S2 dot
    dot(E2) == beta22.*S2.*(I2./N2) +  beta21.*S2.*(I1./N1) - sigma.*E2 ; %E2 dot
    dot(I2) == sigma.*E2 - gamma.*I2 -  phi2.*I2 ; %I2 dot
    dot(R2) == qV.*uV2  + gamma.*I2 - omega.*R2 ; %R2 dot
    dot(N1) == -phi1.*I1 ;
    dot(N2) == -phi2.*I2 ;
    dot(V1) ==  uV1 ;
    dot(V2) ==  uV2 });    
    
    
end

%Constraints on state and control variable
if CASE==1 %%%% No controls

    cbb={icollocate(lowerS <= state)
        icollocate(lowerC == control)
    initial(state == x0ic)};


elseif CASE==3 %%%% Vaccine
    
    cbb={icollocate(lowerS <= state)
        icollocate(0 <= control(1) <= S1 )
        icollocate(0 <= control(2) <= S2)
       icollocate(control(1)+control(2)<=MaxTreat(2).*SCAR)
    initial(state == x0ic)};

elseif CASE==5 %%%% Ad Hoc  Vaccine
    
           % Box constraints 
    cbb={icollocate(lowerS <= state)
        icollocate(0 <= control(1) <= S1 )
        icollocate(0 <= control(2) <= S2 )
        icollocate(control(1)<=(N1./(N1+N2)).*MaxTreat(2).*SCAR)
        icollocate(control(2)<=(N2./(N1+N2)).*MaxTreat(2).*SCAR)
        initial(state == x0ic)};
end  
    

%Objective Function
if OBJ==1 %Cost minimization with Policy Adjustment cost
    
 
 %Policy adjustment cost 
  adj= cAdj.*((N2./(N1+N2)).*uV1 - (N1./(N1+N2)).*uV2).^2 ;   
  

  objective = integrate(exp(-r.*t).*((phi1+w).*cI.*I1 + (phi2+w).*cI.*I2 + cV.*(uV1+uV2) + adj ));
 
end

if OBJ==2
    
     %Policy adjustment cost 
  adj= cAdj.*((N2./(N1+N2)).*uV1 - (N1./(N1+N2)).*uV2).^2 ;   
  

  objective = integrate(exp(-r.*t).*((phi1+w).*cI.*I1 + cV.*(uV1) + adj ));
    
    
end



% Giving a name to each problem to the solver
options=struct;

 if CASE==1
        options.name='No Controls';
    elseif CASE==3
        options.name='Optimal Constrained Vaccine';
    elseif CASE==5
        options.name='Ad Hoc Constrained Vaccine';
end


% Choosing Solver
   options.solver = 'knitro'; % choose solver ('knitro' is other good one)



% Initial Guess
if isempty(GUESS)==1
        %Creating initial guess: Initial guess are initial conditions
    x0_guess = { icollocate({S1 == x0ic(1);
    E1 == x0ic(2);
    I1 == x0ic(3);
    R1 == x0ic(4);
    S2 == x0ic(5);
    E2 == x0ic(6);
    I2 == x0ic(7);
    R2 == x0ic(8);
    N1 == x0ic(9);
    N2 == x0ic(10);
    V1 == x0ic(11);
    V2 == x0ic(12);}) ;
    collocate(uV1== x0ic(1)*100) %Here I assume the initial guess of the controls will be proportional to the amount of S1, and S2
    collocate(uV2== x0ic(5)*100)};


elseif isempty(GUESS)==0 %isempty(GUESS)~=1 
        

    x0_guess = { icollocate({S1 == GUESS(CASE).S1;
        E1==GUESS(CASE).E1;
        I1==GUESS(CASE).I1;
        R1==GUESS(CASE).R1;
        N1==GUESS(CASE).N1;
        S2==GUESS(CASE).S2;
        E2==GUESS(CASE).E2;      
        I2==GUESS(CASE).I2;
        R2==GUESS(CASE).R2;
        N2==GUESS(CASE).N2;
        V1==GUESS(CASE).V1;
        V2==GUESS(CASE).V2;}) ;
        collocate(uV1==GUESS(CASE).uV1)
        collocate(uV2==GUESS(CASE).uV2)};
        
end


  % Solving the problem
  [solution, result]= ezsolve(objective, {ode,cbb,cterm},x0_guess, options);
 
  % This gives the solver one more try using the results from the prior run. this is robust to do before entering your loop
  [solution, result]= ezsolve(objective, {ode,cbb,cterm},solution, options);
        
     
% Changing solver if exitflag~=0
 if true
counter=0;
    while result.ExitFlag~=0 && counter<1  % limits the lopp size

        options.solver = 'snopt';
        [solution, result]= ezsolve(objective, {ode,cbb,cterm},x0_guess, options);
        
        if result.ExitFlag~=0
            
        options.solver = 'npsol';
        [solution, result]= ezsolve(objective, {ode,cbb,cterm},x0_guess, options);
        
        end
        
        if result.ExitFlag~=0
            
        options.solver = 'snopt';
        [solution, result]= ezsolve(objective, {ode,cbb,cterm},x0_guess, options);

        end
        
        counter=counter+1;
    end
    
 end      
   
    
% Storing the results if solver found a good solution
if result.ExitFlag==0
        Results=result; %%% not sure about this
        ts=subs(icollocate(t), solution);
        S1s=subs(icollocate(S1), solution);
        E1s=subs(icollocate(E1), solution);
        I1s=subs(icollocate(I1), solution);
        R1s=subs(icollocate(R1), solution);
        N1s=subs(icollocate(N1), solution);
        S2s=subs(icollocate(S2), solution);
        E2s=subs(icollocate(E2), solution);
        I2s=subs(icollocate(I2), solution);
        R2s=subs(icollocate(R2), solution);
        N2s=subs(icollocate(N2), solution);
        V1s=subs(icollocate(V1), solution);
        V2s=subs(icollocate(V2), solution);
        uV1s=subs(icollocate(uV1), solution);   
        uV2s=subs(icollocate(uV2), solution);   

        
    else
        Results=result; 

end