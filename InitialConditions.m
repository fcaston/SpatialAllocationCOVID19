%To simulate initial conditions
clear variables; 
close all;  


[N1,N2,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,r,cD,cV,cAdj,cI,w,params,R0] = covid_Parameters_months()



%x0_1=[0.9999999,.0000001,0,0]'; % S(0),E(0),I(0),R(0),N(0)
x0_1=[9999999,1,0,0]'; % S(0),E(0),I(0),R(0),N(0)
x0=[x0_1*N1];

N=x0_1(1)+x0_1(2);
    

t0=0; tf=8+0.230137*2; %tf=4-0.230137*1;  % for when uncontrolled R0 is 3.15
tspan=linspace(t0,tf,100);%[t0 tf];



[t,y]=ode15s(@(t,x)COVIDeqsInitialConditions(x, beta11, gamma, sigma, phi,omega), tspan, x0,[]);

if true
    
figure
subplot(221)
plot(t,y(:,1)./N);
subplot(222)
plot(t,y(:,2)./N);
subplot(223)
plot(t,y(:,3)./N);
subplot(224)
plot(t,y(:,4)./N);
    
    
end

IC_1=[y(end,1)./N,y(end,2)./N,y(end,3)./N,y(end,4)./N]



t0=0; tf=9-0.2301370*1; %tf=4+0.2301370*0; for when uncontrolled R0 is 3.15
tspan=linspace(t0,tf,100);%[t0 tf];

[t,y]=ode15s(@(t,x)COVIDeqsInitialConditions(x, beta11, gamma, sigma, phi,omega), tspan, x0,[]);


if true
    
figure
subplot(221)
plot(t,y(:,1)./N);
subplot(222)
plot(t,y(:,2)./N);
subplot(223)
plot(t,y(:,3)./N);
subplot(224)
plot(t,y(:,4)./N);
    
    
end


IC_2=[y(end,1)./N,y(end,2)./N,y(end,3)./N,y(end,4)./N]









        
        function dx=COVIDeqsInitialConditions(x, beta11, gamma, sigma, phi, omega)
        dx = zeros(4,1);
        n1 = x(1)+x(2)+x(3)+x(4); 

        dx(1) = omega.*x(4) - beta11.*x(1).*(x(3)./n1); %S1  - beta12.*x(1).*(x(7)./n2)
        dx(2) = beta11.*x(1).*(x(3)./n1)  - sigma.*x(2)  ; %E1 dot + beta12.*x(1).*(x(7)./n2)       
        dx(3) = sigma.*x(2) - gamma.*x(3) - phi.*x(3)  ; %I1 dot %- phi.*x(3) 
        dx(4) = gamma.*x(3) - omega.*x(4) ; %R1 dot    
 
        end
    
     