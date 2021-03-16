function [N1,N2,omega,beta11,beta22,beta12,beta21,sigma,gamma,phi1,phi2,qV,r,cV,cAdj,cI,w,params,R0] = covid_Parameters_months()

% Population size
params.N1=1e+0;
params.N2=1e+0;


%Epidemiological/Biological parameters


params.sigma = 1/(3/30.4167); %3 day latency, converted to months (Davies et al., 2020)
params.gamma = 1/(5/30.4167);  %5 day recovery, converted to months (Davies et al., 2020)
params.omega=1/6; % 6-month immunity (Edridge et al., 2020)
params.w= 0.133; %Disability weight associated with covid-19 (Nurchis et al., 2020)


params.phi1 =params.gamma./(1/0.0278 -1); % Based on a 0.0178 case fatality rate (Abdollahi et al., 2020), using the fact that phi/(phi+gamma) = 0.0178
params.phi2 =params.gamma./(1/0.0178 -1); % Based on a 0.0178 case fatality rate (Abdollahi et al., 2020), using the fact that phi/(phi+gamma) = 0.0178



params.beta11 = (params.gamma+params.phi1)*2.2*(0.65) ; %Li et al. (2020) for R0 estimate 
params.beta22 = (params.gamma+params.phi2)*2.2*(0.65); %Tian et al. (2020) for nonphamaceutical intervention effect on R0
params.beta12 = (params.gamma+params.phi1)*2.2*(0.97-0.65);
params.beta21 = (params.gamma+params.phi2)*2.2*(0.97-0.65);


%%% Efficiency paramters
params.qV = 0.65; %Efficiency of Vaccines


%Economic parameters
params.r=.015./12; %Discount rate
params.cV = 20*2; % Cost of treating one individual via Vaccine
params.cAdj= 10000;% Policy Adjustment Cost


%%% Value of Statistical Life (VSL)
params.cI= 10e6; 


%Naming the parameters
N1=params.N1; N2=params.N2; 
beta11=params.beta11; beta22=params.beta22; beta12=params.beta12; beta21=params.beta21;
sigma=params.sigma; gamma=params.gamma; omega=params.omega; phi1=params.phi1; phi2=params.phi2;
qV=params.qV;
r=params.r; cV=params.cV ; cAdj=params.cAdj; 
w=params.w; 
cI=params.cI ; 


%Calculating the Next-Generation Matrix and the R0
F=  [beta11 0 beta12 0; 0 0 0 0; beta21 0 beta22 0; 0 0 0 0];
V= [0 sigma 0 0; -(gamma+phi1) -sigma 0 0; 0 0 0 sigma; 0 0 -(gamma+phi2) -sigma];

FV=-F*inv(V);
R0=eigs(FV,1,'lr')
R0_all=eig(FV);


%%% References for Epidemiological Parameters
%Abdollahi, E., Champredon, D., Langley, J. M., Galvani, A. P., and Moghadas, S. M. (2020). Temporal estimates of case-fatality rate for COVID-19 outbreaks in Canada and the United States. CMAJ.
%Davies, N. G., Klepac, P., Liu, Y., et al. (2020). Age-dependent effects in the transmission and control of COVID-19 epidemics. Nature Medicine, 26(8):1205–1211.
%Edridge, A. W., Kaczorowska, J. M., Hoste, A. C., Bakker, M., Klein, M., Jebbink, M. F., Matser, A., Kinsella, C., Rueda, P., Prins, M., et al. (2020). Coronavirus protective immunity is short-lasting. MedRxiv.
%Li, Q., Guan, X., Wu, P., Wang, X., Zhou, L., Tong, Y., Ren, R., Leung, K. S., Lau, E. H., Wong, J. Y., et al. (2020). Early transmission dynamics in wuhan, china, of novel coronavirus–infected pneumonia. New England Journal of Medicine.
%Nurchis, M. C., Pascucci, D., Sapienza, M., Villani, L., D’Ambrosio, F., Castrini, F., Specchia, M. L., Laurenti, P., and Damiani, G. (2020). Impact of the burden of COVID-19 in italy: Results of disability-adjusted life years (dalys) and productivity loss. International Journal of Environmental Research and Public Health, 17(12):4233.
%Tian, H., Liu, Y., Li, Y., Wu, C.-H., Chen, B., Kraemer, M. U., Li, B., Cai, J., Xu, B., Yang, Q., et al. (2020). An investigation of transmission control measures during the first 50 days of the COVID-19 epidemic in china. Science, 368(6491):638–642.

end
