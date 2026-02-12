%% Clear workspace
close all; clear all;

%% Load daTsoilset
addpath("C:\github\GHG-data-useful-scripts\Matlab-scripts")

VE_ori = readtable('sample2.txt'); % ##### Input ##### Please change the file name accordingly
colors = lines(1); % Colors

%% =========================
% 1. Night-time NEE vs Soil Temperature
%% =========================
figure('Position',[100 100 800 800])

hold on
Tsoil = VE_ori.T10; %Either Tsoil or Ta. Choose the one that best fit to the curve
NEE = VE_ori.NEE ;
NEE(NEE==-9999) = NaN;
NEE(VE_ori.PAR > 0.5) = NaN; % night-time only

% Remove NaN and non-positive NEE
idx = ~isnan(Tsoil) & ~isnan(NEE) & NEE>0;
Tsoil = Tsoil(idx);
NEE = NEE(idx);

% Sort
[TsoilSort, sortIdx] = sort(Tsoil);
NEESort = NEE(sortIdx);

% Non-linear model: simple exponential
opts = statset('Display','iter','TolFun',1e-10);
modelFun = @(b,x)(b(1).*exp(b(2).*(1./56.02-1./(x(:,1)+273.15-227.13))));
beTsoil0 = [1.2,337];
mdl = fitnlm(TsoilSort, NEESort, modelFun, beTsoil0,'Options',opts);

% Prediction
TsoilFit = linspace(min(TsoilSort), max(TsoilSort), 100)';
NEE_fit = predict(mdl, TsoilFit);

% Plot
plot(TsoilFit, NEE_fit, 'LineWidth',2,'Color',colors);
scatter(TsoilSort, NEESort, 10, colors, 'filled', 'MarkerFaceAlpha',0.2);

% ===== Show parameters in command window =====
disp('Night-time Reco fit parameters:')
disp(mdl.Coefficients)
R2_reco = mdl.Rsquared.Ordinary;
disp(['Reco fit R2 = ' num2str(R2_reco)])

xlabel('Soil Temperature (^oC)');
ylabel('Night-time NEE (\mumol CO_2 m^{-2} s^{-1})');
box on; ax=gca; ax.FontSize=14; ax.LineWidth=1;



hold off
%% =========================
% 2. PAR vs GPP
%% =========================
figure('Position',[100 100 800 800])
hold on
PAR = VE_ori.PAR ; 
NEE = VE_ori.NEE;
NEE(NEE==-9999) = NaN;
% Estimate GPP: GPP = Reco - NEE (simplified)
% Fit Reco first (night-time)
NEE_night = NEE;
NEE_night(VE_ori.PAR > 0.5) = NaN;
idxNight = ~isnan(NEE_night);
TsoilNight = VE_ori.T10(idxNight);
NEE_night = NEE_night(idxNight);
opts = statset('Display','iter','TolFun',1e-10);
modelFunReco = @(b,x)(b(1).*exp(b(2).*(1./56.02-1./(x(:,1)+273.15-227.13))));
beTsoil0 = [1.2,337];
mdlReco = fitnlm(TsoilNight, NEE_night, modelFunReco, beTsoil0,'Options',opts);

% Predict Reco for all PAR points
Reco_est = predict(mdlReco, VE_ori.T10);
GPP_est = Reco_est - NEE;

% Remove invalid
idx = ~isnan(PAR) & ~isnan(GPP_est);
PAR = PAR(idx);
GPP_est = GPP_est(idx);

% Fit Michaelis-Menten
mmFun = @(b,x) b(1).*x.*b(2)./(b(1).*x + b(2));
beTsoil0 = [-0.1 -25];
mdlGPP = fitnlm(PAR, GPP_est, mmFun, beTsoil0);

% Prediction
PAR_fit = linspace(min(PAR), max(PAR), 100)';
GPP_fit = predict(mdlGPP, PAR_fit);

% Plot
plot(PAR_fit, GPP_fit, 'LineWidth',2,'Color',colors);
scatter(PAR, GPP_est, 10, colors, 'filled', 'MarkerFaceAlpha',0.2);

% ===== Show parameters in command window =====
disp('GPP light-response fit parameters:')
disp(mdlGPP.Coefficients)
R2_gpp = mdlGPP.Rsquared.Ordinary;
disp(['GPP fit R2 = ' num2str(R2_gpp)])

xlabel('PAR (\mu mol m^{-2} s^{-1})');
ylabel('GPP (\mumol CO_2 m^{-2} s^{-1})');
box on; ax=gca; ax.FontSize=14; ax.LineWidth=1;
hold off