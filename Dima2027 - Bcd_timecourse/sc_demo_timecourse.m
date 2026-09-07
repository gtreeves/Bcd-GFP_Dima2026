% sc_demo_timecourse
%
% Run this script for the demonstration of the time course RICS code
% (tRICS) on a single czi acquisition, not using tiles.
%
% This script is distributed on the CC BY-SA license, created by Greg
% Reeves 2026, Chemical Engineering Department, Texas A&M University

clear
close all

%
% Basic input parameters
%
filename = 'demo_czi_files\Timecourse\nc14-1.czi';
data_ch = 1; % data channel, in which Bcd-GFP is present
genotype = 'Bcd-GFP';
sws = 5;
mask_ch = 2; % mask channel, in which His2Av-RFP is present
usemask = true;
istart = []; % index of the czi time course file that we start with
iend = []; % index of the czi time course file that we end with
basename = 'embryo01'; % if the output file needs a base name
zerospad = false; % option to pad data with zeros to prevent ACF wraparound
hasbackground = false; % if some of the image frame is outside the embryo
yesplot = true; % if true, outputs diagnostic images


%
% Geometric params
%
w0 = 0.25;
wz = 0.75;
chi = wz/w0;
w0r = 0.3; % ad hoc
gamma = sqrt(2)/4;
V_eff = pi^1.5*chi*w0^3/gamma;
Vnuc_ch = pi^1.5*chi*w0r^3/gamma;
w0_cc = sqrt(0.5*(w0.^2 + w0r.^2));
Vcc = pi^1.5*chi*w0_cc.^3/gamma;



%
% Calculate timecourse of ACFs
%
[~,Gs_t] = analyze_RICS_t(filename,data_ch,genotype,w0,wz,...
			sws,mask_ch,usemask,istart,iend,basename,...
			zerospad,hasbackground,yesplot);

%
% Fitting model to ACFs
%
B0 = [0 -1e-3 1e-3]; % set B0 to zero so we don't overfit if it drops below zero
D20 = 3; % choose plausible value of the unidentifiable diffusivity
W0 = []; % allow w0 to vary slightly for robustness of fit
data = fit_both_models_2step(Gs_t,B0,D20,W0);



%
% Save calcultions to mat files
%
save(['Mat/',outfilename,'_timecourse_Gs_t'],'Gs_t','data','-v7.3') % ACFs (large file)
save(['Mat/',outfilename,'_timecourse_data'],'data') % tertiary data




%
% Calculate concentrations and fractions from the estimated parameters
%
t = data.t;
Cnuc = 1./(V_eff*data.Anuc)/ 0.602; % converted to nM
phi = data.phinuc;
phicorrel = Vcc./Vnuc_ch.*data.Acc'./data.Anuc_ch;
phiuncorrel = phi - phicorrel;
phifree = 1 - phi;
Cnuc_ch = 1./(V_eff*data.Anuc_ch)/ 0.602; % converted to nM

%
% Make plots
%
figure('pos',[9         102        1899         592])
TL = tiledlayout(2,5,'TileSpacing','Compact','Padding','Compact');

%
% First row
%
nexttile
plot(t,Cnuc)
xlabel('time [min]')
ylabel('tot conc [nM]')

nexttile
plot(t,Cnuc.*phifree)
xlabel('time [min]')
ylabel('free conc [nM]')

nexttile
plot(t,Cnuc.*phicorrel)
xlabel('time [min]')
ylabel('correl conc [nM]')

nexttile
plot(t,Cnuc.*phiuncorrel)
xlabel('time [min]')
ylabel('uncorrel conc [nM]')

nexttile
plot(Cnuc.*phifree,Cnuc.*phicorrel)
xlabel('free conc [nM]')
ylabel('correl conc [nM]')

%
% second row
%
nexttile
plot(t,Cnuc_ch)
xlabel('time [min]')
ylabel('H2Av conc [nM]')

nexttile
plot(t,phifree)
xlabel('time [min]')
ylabel('free frac [nM]')

nexttile
plot(t,phicorrel)
xlabel('time [min]')
ylabel('correl frac [nM]')

nexttile
plot(t,phiuncorrel)
xlabel('time [min]')
ylabel('uncorrel frac [nM]')

nexttile
plot(Cnuc.*phifree,Cnuc.*phiuncorrel)
xlabel('free conc [nM]')
ylabel('uncorrel conc [nM]')










