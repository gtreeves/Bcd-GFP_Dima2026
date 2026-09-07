% sc_demo_rectangle
%
% This script analyze one Bcd-GFP file that is in a rectangular shape
% (4096x1024, or the transpose). It will divide the rectangle into four
% squares and calculate ACFs (secondary data) and estimate parameters
% (tertiary data) for each square.
%
% This script is distributed on the CC BY-SA license, created by Greg
% Reeves 2026, Chemical Engineering Department, Texas A&M University

clear
close all

%
% Basic input parameters
%
filename = 'demo_czi_files\Rectangle\nc14-1-l.czi'; % rectangle with long line time
% alternate demo:
% filename = 'demo_czi_files\Rectangle\nc14-1-s.czi'; % rectangle with short line time
data_ch = 1;
genotype = 'Bcd-GFP';
subtravg = true;
sws = 5;
mask_ch = 2;
usemask = 12;
fitprofiles = true;
outfilename = 'Bcd_rectangle';
zerospad = false;
hasbackground = false;
yesplot = true;


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
% Location data (specific to this embryo)
%
d = 180; % location of left side of tile in microns from anterior pole
L = 504; % embryo length
orientation = 1; % anterior to the left

% calculation of locations of centers of each square
dr = 0.0325; % pixel size (all images taken at dr = 0.0325 micron/pixel)
S0 = (d + dr*(1024*[0 1 2 3] + 512)).*(orientation == 1) + ...
	(d + dr*(1024*[3 2 1 0] + 512)).*(orientation == -1) + ...
	d.*isnan(orientation); % 19-by-4

%
% Metadata
%
w0 = 0.25; % Used to be 0.24
wz = 0.75;


istart = [];
iend = [];
xilimit = 51;

%
% Calculate ACFs
%
Gs = analyze_RICS_x(filename,data_ch,genotype,w0,wz,...
	sws,mask_ch,usemask,zerospad,hasbackground,xilimit,yesplot);

%
% Fitting model to ACFs
%
B0 = [0 -1e-3 1e-3]; % set B0 to zero so we don't overfit if it drops below zero
D20 = 3; % just for the 2-cpt model for now
W0 = []; % allow w0 to vary slightly for robustness of fit
data = fit_both_models_2step(Gs,B0,D20,W0);

%
% Save calcultions to mat files
%
save(['Mat/',outfilename,'_rectangle_Gs'],'Gs','data','-v7.3') % ACFs (large file)
save(['Mat/',outfilename,'_rectangle_data'],'data') % tertiary data


%
% Make plots
%
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
C = colormap('lines');

%
% First row
%
nexttile
plot(S0,Cnuc)
xlabel('AP coord [\mum]')
ylabel('tot conc [nM]')

nexttile
plot(S0,Cnuc.*phifree)
xlabel('AP coord [\mum]')
ylabel('free conc [nM]')

nexttile
plot(S0,Cnuc.*phicorrel)
xlabel('AP coord [\mum]')
ylabel('correl conc [nM]')

nexttile
plot(S0,Cnuc.*phiuncorrel)
xlabel('AP coord [\mum]')
ylabel('uncorrel conc [nM]')

nexttile
plot(Cnuc.*phifree,Cnuc.*phicorrel)
xlabel('free conc [nM]')
ylabel('correl conc [nM]')

%
% second row
%
nexttile
plot(S0,Cnuc_ch)
xlabel('AP coord [\mum]')
ylabel('H2Av conc [nM]')

nexttile
plot(S0,phifree)
xlabel('AP coord [\mum]')
ylabel('free frac [nM]')

nexttile
plot(S0,phicorrel)
xlabel('AP coord [\mum]')
ylabel('correl frac [nM]')

nexttile
plot(S0,phiuncorrel)
xlabel('AP coord [\mum]')
ylabel('uncorrel frac [nM]')

nexttile
plot(Cnuc.*phifree,Cnuc.*phiuncorrel)
xlabel('free conc [nM]')
ylabel('uncorrel conc [nM]')







