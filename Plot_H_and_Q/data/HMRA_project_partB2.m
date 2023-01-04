% F:\POLITO\Hydro-meteorological risk assessment\Project
clear all ;
clc; 
format short g, 
close all;

%% read data:
T = readtable('date.txt','DatetimeType','datetime');
date = datetime(table2array(T));
Q = importdata("Q.txt");

%A1
%% Compute the annual maximum discharge:
test = double.empty;
Qmaxima_year = double.empty;
for j=1951:2019
    ind = find(year(date) == j);
    test = Q(ind);
    maxima = max(test);
    Qmaxima_year = [Qmaxima_year ; maxima];
end
%% plot time series
years= 1951:2019;
figure(1)
plot(years, Qmaxima_year,'-*')
xlabel("Years")
ylabel("Maximum annual Discharge [m³/s]")
%% sample description
Q_mean=mean(Qmaxima_year); %mu
Q_var=var(Qmaxima_year); %sigma squared, s²
Q_std=std(Qmaxima_year); %sigma, s
Q_moment= moment(Qmaxima_year,3);
Q_max=max(Qmaxima_year);
CV = Q_std/Q_mean;
M3 = moment(Q_max, 3);
CA = M3/Q_std^3;
%% Frequency histogram
k=round(2*length(Qmaxima_year)^0.4);
figure (2)
histogram(Qmaxima_year,k,'Normalization','Probability')
title ('Frequency Histogram')
xlabel ('Discharge Qmax [m^3/s]')
ylabel ('k');
%% empirical non-exceedance frequency function
Qmax_sorted=sort(Qmaxima_year);
n = length(Qmax_sorted);
Fi= ((1:n)/(n+1));
figure(3), ecdf(Qmax_sorted)
title('Empirical Non-Exceedance Frequency Function');

%A2
%% exponential probability distribution
Qmax_sorted = sort(Qmaxima_year, 'ascend');
Fi = (1:n) ./ (n+1) ;
ui_exp = -log(1-Fi);
figure(4), plot(Qmax_sorted,ui_exp,'*')
title('Probability plot - Exponential Distribution')
xlabel('Discharge Qmax [m^3/s]'), ylabel('Reduced Variate [-]')

%% Gumbel distribution
Qmax_sorted = sort(Qmaxima_year, 'ascend');
Fi = (1:n) ./ (n+1) ;
ui_gumb = -log(-log(Fi));
figure(5), plot(Qmax_sorted,ui_gumb,'*')
title('Probability plot - Gumbel Distribution')
xlabel('Discharge Qmax [m^3/s]'), ylabel('Reduced Variate [-]')
%% normal distribution
figure(6), plot(Qmax_sorted,norminv(Fi),'*')
title('Probability plot - Normal Distribution')
xlabel('Discharge Qmax [m^3/s]'), ylabel('Reduced Variate [-]')

%% lognormal distribution
yi=log(Qmax_sorted);
figure(7), plot(yi,norminv(Fi),'*')
title('Probability plot - LogNormal Distribution')
xlabel('Discharge Qmax [m^3/s]'), ylabel('Reduced Variate [-]')

%% GEV distribution
ui_gev = -log(Fi);%reduced variate, is it right like this?
figure(), plot(Qmax_sorted,ui_gev,'*')
title('Probability plot - GEV Distribution')
xlabel('Discharge Qmax [m^3/s]'), ylabel('Reduced Variate [-]')


%% A3 
%{
Estimate the parameters for those distributions that successfully comply 
with the preliminary assessments plus the GEV distribution. Estimate the 
parameters with the methods of moments and L-moments. Where possible, 
plot the distributions with the estimated parameters into the probability
plots
%}

%% estimation of parameters for exponential distribution, 
% method of moments (mom)
theta_exp = Q_mean; 
ui_expMom = Qmax_sorted/theta_exp;

figure(8), plot(Qmax_sorted,ui_exp,'*')
title('Exponential Probability Plot with Estimated Parameter by Method of Moments')
xlabel('Xi'), ylabel('Reduced Variate [-]')
hold on
plot(Qmax_sorted,ui_expMom,'r-')
hold off
%{
%method of L-moments: 
L1 = theta_exp; 
L2 = theta_exp/2;
tau3 = 1/3;
tau4 = 1/6; 
%}


%% estimation of parameters for Gumbel distribution, mom
theta1_gumMom = Q_mean-0.5772*Q_std*(sqrt(6)/pi);
theta2_gumMom = Q_std*(sqrt(6)/pi);
ui_gumMom = (Qmax_sorted-theta1_gumMom)/(theta2_gumMom);
figure(8), plot(Qmax_sorted,ui_gumb,'*')
title('Gumbel Probability Plot with Estimated Parameters by Method of Moments')
%eventually rename this plot?
xlabel('Discharge Qmax [m³/s]'), ylabel('Reduced Variate [-]')
hold on
plot(Qmax_sorted,ui_gumMom,'r-')
hold off

%% estimation of parameters for normal distribution, mom
theta1_norMom = Q_mean;
theta2_norMom = Q_std;
ui_norMom = (Qmax_sorted - theta1_norMom)/theta2_norMom;
figure(9), plot(Qmax_sorted,norminv(Fi),'*')
title('Normal Probability Plot wit-h Estimated Parameters by Method of Moments')
xlabel('Discharge Qmax [m³/s]'), ylabel('Reduced Variate [-]')
hold on
plot(Qmax_sorted,ui_norMom,'r-')
hold off


%% estimation of parameters for lognormal distribution, mom
theta1_logMom = log(Q_mean) - 1/2 * log(1 + Q_var/(Q_mean .^2));
theta2_logMom = sqrt(log(1 + Q_var/(Q_mean .^2)));
ui_logMom = (yi -theta1_logMom)/theta2_logMom;
figure(10), plot(yi, norminv(Fi), '*')
title('Lognormal Probability Plot with Estimated Parameters by Method of Moments')
xlabel('Discharge Qmax [m³/s]'), ylabel('Reduced Variate [-]')
hold on 
plot(yi, ui_logMom, 'r-')
hold off

%Copied code from other report and adjusted to our variables
logp = log(Qmax_sorted); % = yi
teta1lgmethodofmoments = mean(logp);
teta2lgmethodofmoments = sqrt(var(logp));
figure(11), plot(logp,norminv(Fi),'*')
title('LogNormal Probability plot with estimated parameters by method of moments')
xlabel('Xi'), ylabel('Reduced Variate(Ui)')
ui_analyticallognormalmoment = (logp-teta1lgmethodofmoments)/(teta2lgmethodofmoments);
hold on
plot(logp,ui_analyticallognormalmoment,'r-')
hold off

%slightly different results, theta 1 is almost the same, theta 2 is a bit
%different, data aligns a bit differently

%% estimation of parameters for GEV distribution, mom
tt3 = 0.01:0.001:0.03; %what is this and where does it come from? 
% -> initial estimation of theta3?

%calculation of CA with theta3 and the gamma-function
cavgev = abs(tt3)./tt3 .* (-gamma(1+3*tt3)+3*gamma(1+tt3).*gamma(1+2*tt3) -2*gamma(1+tt3).^3)./ (gamma(1+2*tt3)-gamma(1+tt3).^2).^(3/2);
y = cavgev-CA; %difference between the estimated and the observed CA
figure(12) ,plot(tt3,y)
xlabel('Ɵ3')
ylabel('\Delta CA')
title('Graphical Derivation of \Theta3')
hold on, plot (tt3,tt3*0,'k-'), hold off

tt3gevok= 0.019;
tt2gevok= sqrt((Q_var*tt3gevok^2)/((gamma(1+2.*tt3gevok))-((gamma(1+tt3gevok)^2))));
tt1gevok=Q_mean-((tt2gevok*(1-gamma(1+tt3gevok)))/tt3gevok);

% ui_gevMom = exp(-(1-(tt3gevok/tt2gevok*(Qmax_sorted-tt1gevok))).^(1/tt3gevok));
% figure(13), plot(Qmax_sorted, ui_gev,'*')
% title('GEV Probability plot with estimated parameters by method of moments')
% xlabel('Discharge Qmax [m^3/s]'), ylabel('Reduced variate [-]')
% hold on
% plot(Qmax_sorted,ui_gevMom,'r-')
% hold off
%% not right like this

%no idea what is going on here or how it is supposed to help us

%% Method of L-Moments
%general
b0 = Q_mean;
nQ = length(Qmaxima_year);

i = 1:nQ; 
w = (i-1)/(nQ-1);
arg = w .* Qmax_sorted; 

b1 = 1/nQ * sum((i-1)/(nQ-1)*Qmax_sorted);
b2 = 1/nQ * sum(((i-1).*(i-2))/((nQ-1).*(nQ-2))*Qmax_sorted);

l1 = Q_mean;
l2 = 2*b1 - b0;
l3 = 6*b2-6*b1+b0;
%{
%a second waay to compute the above:
u= 0 ;
for k = 1 : nQ
    u = u + ((k-1)/(nQ-1))*Qmax_sorted(k,1);
end 
b1 = u/nQ;

uu=0;
for kk = 1 : nQ
    uu = uu + (((kk-1)*(kk-2))/((nQ-1)*(nQ-2)))*Qmax_sorted(kk,1);
end 
b2 = uu/nQ
%}

%% parameters for exponential distribution Method of L-Moments (Lm)
theta_expLm = l1;
ui_expLm = Qmax_sorted/theta_expLm;
figure(13), plot(Qmax_sorted,ui_exp,'*')
title('Exponential Probability Plot with Estimated Parameter by Method of L-Moments')
xlabel('DIscharge Qmax [m³/s]'), ylabel('Reduced Variate [-]')
hold on
plot(Qmax_sorted,ui_expLm,'r-')
hold off

%% parameters for Gumbel distribution Method of L-Moments (Lm)
theta2_gumLm = l2/log(2);
theta1_gumLm = l1 - 0.5772 * theta2_gumLm;
ui_gumLm = (Qmax_sorted - theta1_gumLm)/(theta2_gumLm);
figure(14), plot(Qmax_sorted,ui_gumb,'*')
title('Gumbel Probability Plot with Estimated Parameters by Method of L-Moments')
xlabel('Discharge Qmax [m³/s]'), ylabel('Reduced Variate [-]')
hold on
plot(Qmax_sorted,ui_gumLm,'r-')
hold off

%% parameters for Normal distribution Method of L-Moments (Lm)
theta1_norLm = l1;
theta2_norLm = power(pi, 0.5) * l2;
ui_norLm = (Qmax_sorted - theta1_norLm)/theta2_norLm;
figure(15), plot(Qmax_sorted,norminv(Fi),'*')
title('Normal Probability Plot with Estimated Parameters by Method of L-Moments')
xlabel('Discharge Qmax [m³/s]'), ylabel('Reduced Variate [-]')
hold on
plot(Qmax_sorted,ui_norLm,'r-')
hold off

%% parameters for Lognormal distribution Method of L-Moments (Lm)
theta2_logLm = sqrt(2)*norminv((1+l2/l1)/2);
theta1_logLm = log(l1) - power(theta2_logLm, 2)/2;
ui_logLm = (yi -theta1_logLm)/theta2_logLm;
figure(16), plot(yi,norminv(Fi),'*')
title('Lognormal Probability Plot with Estimated Parameters by Method of L-Moments')
xlabel('Discharge Qmax [m³/s]'), ylabel('Reduced Variate [-]')
hold on 
plot(yi, ui_logLm, 'r-')
hold off

%% parameters for GEV distribution Method of L-Moments (Lm)
c = 2/(3+l3/l2) - log(2)/log(3);
theta3_gevLm = 7.859*c + 2.9554*c^2;
theta2_gevLm = l2*theta3_gevLm / ((1-power(2, -theta3_gevLm))* gamma(1+theta3_gevLm));
theta1_gevLm = l1 - theta2_gevLm/theta3_gevLm * (1-gamma(1+theta3_gevLm));

%do we have to do something else with the GEV parameters??

%B1
%% Pearson test parameter
x = Qmax_sorted;
n = numel(x); %sample size
k = round(2*(n^0.4)); %number of classes
alpha = 0.05; %significance level
Ei = n * 1/k; %expected frequecies in the i-th class
P_x = (0:k)./k; %probability


%% Pearson test - Lognormal, L-Moments
x_P_logLm = (norminv(P_x)*theta2_logLm+theta1_logLm);
for classes = 1:k
xi_logLm = logical((x >= x_P_logLm(classes)) & (x < x_P_logLm(classes+1)));
Oi_logLm(classes) = sum(xi_logLm);
end
X2_logLm = sum (((Oi_logLm(1:k)-Ei).^2)./Ei);
X2lim_logLm= chi2inv(1-alpha,k-2-1);

if X2_logLm <= X2lim_logLm

    disp ('Test passed, H0 accepted, Lognormal distribution parameters by method of L-moments has good fitting')
else
    disp ('Test NOT passed, H0 rejected, Lognormal distribution parameters by method of L-moments is not suitable')
end

%% Pearson test - Lognormal, Moments
x_P_logMom = (norminv(P_x)*theta2_logLm+theta1_logLm);
for classes = 1:k
xi_logMom = logical((x >= x_P_logMom(classes)) & (x < x_P_logMom(classes+1)));
Oi_logMom(classes) = sum(xi_logMom);
end
X2_logMom = sum (((Oi_logMom(1:k)-Ei).^2)./Ei);
X2lim_logMom= chi2inv(1-alpha,k-2-1);

if X2_logMom <= X2lim_logMom

    disp ('Test passed, H0 accepted, Lognormal distribution parameters by method of moments has good fitting')
else
    disp ('Test NOT passed, H0 rejected, Lognormal distribution parameters by method of moments is not suitable')
end

%% Pearson test - Gumble, Moments
x_P_gumMom = log(-log(P_x))*theta2_gumMom+theta1_gumMom
for classes = 1:k
xi_gumMom = logical((x >= x_P_gumMom(classes)) & (x < x_P_gumMom(classes+1)));
Oi_gumMom(classes) = sum(xi_gumMom);
end
X2_gumMom = sum (((Oi_gumMom(1:k)-Ei).^2)./Ei);
X2lim_gumMom= chi2inv(1-alpha,k-2-1);

if X2_gumMom <= X2lim_gumMom

    disp ('Test passed, H0 accepted, Gumble distribution parameters by method of moments has good fitting')
else
    disp ('Test NOT passed, H0 rejected, Gumble distribution parameters by method of moments is not suitable')
end

%% Pearson test - Gumble, L-Moments
x_P_gumLm = log(-log(P_x))*theta2_gumLm+theta1_gumLm
for classes = 1:k
xi_gumLm = logical((x >= x_P_gumLm(classes)) & (x < x_P_gumLm(classes+1)));
Oi_gumLm(classes) = sum(xi_gumLm);
end
X2_gumLm = sum (((Oi_gumLm(1:k)-Ei).^2)./Ei);
X2lim_gumLm= chi2inv(1-alpha,k-2-1);

if X2_gumLm <= X2lim_gumLm

    disp ('Test passed, H0 accepted, Gumble distribution parameters by method of L-moments has good fitting')
else
    disp ('Test NOT passed, H0 rejected, Gumble distribution parameters by method of L-Moments is not suitable')
end

%% Pearson test - Normal, L-Moments
x_P_norLm = norminv(P_x)*theta2_norLm+theta1_norLm;
for classes = 1:k
xi_norLm = logical((x >= x_P_norLm(classes)) & (x < x_P_norLm(classes+1)));
Oi_norLm(classes) = sum(xi_norLm);
end
X2_norLm = sum (((Oi_norLm(1:k)-Ei).^2)./Ei);
X2lim_norLm= chi2inv(1-alpha,k-2-1);

if X2_norLm <= X2lim_norLm

    disp ('Test passed, H0 accepted, Normal distribution parameters by method of L-Moments has good fitting')
else
    disp ('Test NOT passed, H0 rejected, Normal distribution parameters by method of L-Moments is not suitable')
end

%% Pearson test - Normal, Moments
x_P_norMom = norminv(P_x)*theta2_norMom+theta1_norMom;
for classes = 1:k
xi_norMom = logical((x >= x_P_norMom(classes)) & (x < x_P_norMom(classes+1)));
Oi_norMom(classes) = sum(xi_norMom);
end
X2_norMom = sum (((Oi_norMom(1:k)-Ei).^2)./Ei);
X2lim_norMom= chi2inv(1-alpha,k-2-1);

if X2_norMom <= X2lim_norMom

    disp ('Test passed, H0 accepted, Normal distribution parameters by method of Moments has good fitting')
else
    disp ('Test NOT passed, H0 rejected, Normal distribution parameters by method of Moments is not suitable')
end

%% Pearson test - Exponential, Moments
x_P_expMom = -log(1-P_x)*theta_exp;
for classes = 1:k
xi_expMom = logical((x >= x_P_expMom(classes)) & (x < x_P_expMom(classes+1)));
Oi_expMom(classes) = sum(xi_expMom);
end
X2_expMom = sum (((Oi_expMom(classes)-Ei).^2)./Ei);
X2lim_expMom= chi2inv(1-alpha,k-2-1);

if X2_expMom <= X2lim_expMom

    disp ('Test passed, H0 accepted, Exponential distribution parameters by method of Moments has good fitting')
else
    disp ('Test NOT passed, H0 rejected, Exponential distribution parameters by method of Moments is not suitable')
end

%% Pearson test - Exponential, L-Moments
x_P_expLm = -log(1-P_x)*theta_expLm;
for classes = 1:k
xi_expLm = logical((x >= x_P_expLm(classes)) & (x < x_P_expLm(classes+1)));
Oi_expLm(classes) = sum(xi_expLm);
end
X2_expLm = sum (((Oi_norLm(1:k)-Ei).^2)./Ei);
X2lim_expLm= chi2inv(1-alpha,k-2-1);

if X2_expLm <= X2lim_expLm

    disp ('Test passed, H0 accepted, Exponential distribution parameters by method of L-Moments has good fitting')
else
    disp ('Test NOT passed, H0 rejected, Exponential distribution parameters by method of L-Moments is not suitable')
end


%% B2  Anderson Darling Test

%% Anderson test - Exponential
i = 1:length(Qmax_sorted);
Pexp(i) = 1-exp(-Qmax_sorted/Q_mean);
A2_exp = -n-(1/n)*sum(((2*i-1)*log(Px_exp(i))+(2*n+1-2*i)*log(1-Px_exp(i))));
eta = 1.141;
beta = 0.229;
xi = 0.169;
if 1.2*xi < A2_exp
w_exp = 0.0403 + 0.116*(((A2_exp-xi)/beta)^(eta/0.861));
else
w_exp=(0.0403+0.116*((0.2*xi)/beta).^(eta/0.861)).*((A2_exp-0.2*xi)/xi);
end

%% Anderson test - normal
i = 1:length(Qmax_sorted);
Pnorm(i) = normcdf((Qmax_sorted-Q_mean)/Q_std);
A2_norm = -n-(1/n)*sum(((2*i-1).*log(Pnorm(i))+(2*n+1-2*i).*log(1-Pnorm(i))));
eta = 1.147;
beta = 0.229;
xi = 0.167;
if 1.2*xi < A2_norm
w_norm = 0.0403 + 0.116*(((A2_norm-xi)/beta)^(eta/0.861));
else
w_norm=(0.0403+0.116*((0.2*xi)/beta).^(eta/0.861)).*((A2norm-0.2*xi)/xi);
end







