

clear variables; 
%close all;  


[N1,N2,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi,qD,qV,r,cD,cV,cAdj,cI,w,params,R0] = covid_Parameters_months()



x0_1=[0.999999,.000001,0,0,0.999999,.000001,0,0]'; % S(0),E(0),I(0),R(0),N(0)

x0=[x0_1*N1];



time=12;
t0=0; tf=1*time; 
tspan=linspace(t0,tf,100);%[t0 tf];




[t,y]=ode15s(@(t,x)COVIDeqsInitialConditions(x, beta11, beta12, beta22, beta21, gamma, sigma, phi,omega), tspan, x0,[]);

if true
    
figure
subplot(221)
plot(t,y(:,1));
subplot(222)
plot(t,y(:,2));
subplot(223)
plot(t,y(:,3));
subplot(224)
plot(t,y(:,4));
    
    
end


Constraint_1=max(y(:,3))


beta11=beta22;
beta12=0;

[t,y]=ode15s(@(t,x)COVIDeqsInitialConditions(x, beta11, beta12, beta22, beta21, gamma, sigma, phi,omega), tspan, x0,[]);

if true
    
figure
subplot(221)
plot(t,y(:,1));
subplot(222)
plot(t,y(:,2));
subplot(223)
plot(t,y(:,3));
subplot(224)
plot(t,y(:,4));
    
    
end


Constraint_2=max(y(:,3))


Constraint=Constraint_1 + Constraint_2



        
        function dx=COVIDeqsInitialConditions(x, beta11, beta12, beta22, beta21, gamma, sigma, phi, omega)
        dx = zeros(8,1);
        n1 = x(1)+x(2)+x(3)+x(4); 
        n2 = x(5)+x(6)+x(7)+x(8); 

        dx(1) = omega.*x(4) - beta11.*x(1).*(x(3)./n1) - beta12.*x(1).*(x(7)./n2); %S1  - beta12.*x(1).*(x(7)./n2)
        dx(2) = beta11.*x(1).*(x(3)./n1) + beta12.*x(1).*(x(3)./n1) - sigma.*x(2)  ; %E1 dot + beta12.*x(1).*(x(7)./n2)       
        dx(3) = sigma.*x(2) - gamma.*x(3) - phi.*x(3)  ; %I1 dot %- phi.*x(3) 
        dx(4) = gamma.*x(3) - omega.*x(4) ; %R1 dot    
        
        dx(5) = omega.*x(8) - beta22.*x(5).*(x(7)./n2) - beta21.*x(5).*(x(3)./n1); %S1  - beta12.*x(1).*(x(7)./n2)
        dx(6) = beta22.*x(5).*(x(7)./n2) + beta21.*x(5).*(x(3)./n1) - sigma.*x(6)  ; %E1 dot + beta12.*x(1).*(x(7)./n2)       
        dx(7) = sigma.*x(6) - gamma.*x(7) - phi.*x(7)  ; %I1 dot %- phi.*x(3) 
        dx(8) = gamma.*x(7) - omega.*x(8) ; %R1 dot    
 
        end
    
        
