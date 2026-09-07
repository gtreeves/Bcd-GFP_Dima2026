% sc_demo_AP_dataset
%
% Run this script to load the tertiary data from the time course RICS
% (tRICS) analysis during nc14 at different AP coordinates, then plot the
% results. 
%
% This script is distributed on the CC BY-SA license, created by Greg
% Reeves 2026, Chemical Engineering Department, Texas A&M University


clear
close all

%
% Load timecourses
%
load Mat\2023-10-13_20-29-53_Bcd_timecourse Soln
Soln_orig = Soln;
load Mat\2025-05-02_16-20-06_anterior
SolnA = Soln;
load Mat\2025-05-02_16-07-08_posterior
SolnP = Soln;

%
% Geometric params
%
w00 = 0.25;
w0r = 0.3;
chi = 4;
gamma = sqrt(2)/4; % Brown2008 has this as sqrt(2)/4, not sqrt(2)/2
% I am not sure where I got the "sqrt(2)/2" from. So I changed it.
V_eff = pi^1.5*chi*w00^3/gamma;
Vnuc_ch = pi^1.5*chi*w0r^3/gamma;
w0_cc = sqrt(0.5*(w00.^2 + w0r.^2));
Vcc = pi^1.5*chi*w0_cc.^3/gamma;




%% ========================================================================
% Calibrate posterior with anterior
% =========================================================================
% Use the anterior to calibrate the posterior. NOTE: embryos are already
% properly sorted, so SolnA(i) is the same embryo as SolnP(i), so no need
% to go into embryo names 
%
% Calculation of conversion between RICS amplitude and intensity:
% K = A(anterior)*I(anterior)
%
% Then we apply this conversion factor to the posterior: 
% A(posterior) = K/I(posterior)
% {

k = 5; % the fifth column is nc 14
for i = 1:length(SolnA)
	K = mean(SolnA(i).Anuc(:,k).*SolnA(i).nucsignal(:,k),'omitmissing');
	Anuc = SolnP(i).Anuc;
	SolnP(i).Anuc = K./SolnP(i).nucsignal;

end

%
% Concat
%
Soln = [Soln_orig;SolnA;SolnP];
n_embryos = length(Soln);


%}



%% ========================================================================
% Normalize nuclear cycle lengths
%
% Each embryo differs in its exact nuclear cycle length. Here we analyze
% each of embryo and normalize their nc lengths to match each other. See
% Reeves et al., 2012, Dev Cell
%
% =========================================================================
% {


NC_list = 10:14;
tstart = NaN(n_embryos,length(NC_list));
tend = NaN(n_embryos,length(NC_list));
for iii = 1:n_embryos
	t1 = Soln(iii).t;

	t14 = t1(:,end);
	t14_0 = min(t14); % start of nc 14

	if isdatetime(t1)
		tstart(iii,:) = minutes(min(t1) - t14_0);
		tend(iii,:) = minutes(max(t1) - t14_0);
	else
		tstart(iii,:) = min(t1) - t14_0;
		tend(iii,:) = max(t1) - t14_0;
	end

end
tstartmean = mean(tstart,'omitmissing');
tendmean = mean(tend,'omitmissing');

%
% Create general mesh of time points
%
dt = 2; % roughly two minutes between each time point
t_mesh = cell(1,length(NC_list));
Delta_t = tendmean - tstartmean;
nt = round(Delta_t/dt) + 2;
for i = 1:length(NC_list)
	t_mesh{i} = linspace(tstartmean(i),tendmean(i),nt(i))';
end

%}



%% ========================================================================
% Extract AP coordinates of each embryo 
% =========================================================================
% The tabulated AP coordinates of each tRICS acquisition for each embryo.
% {

%
% Distance in microns from anterior pole
%
s = [0.306
	0.280
	0.189
	0.237
	0.164
	0.415
	0.296
	0.425
	0.210
	0.188
	0.126
	0.202
	0.516
	0.178
	0.383
	0.216
	0.339
	0.264
	0.194
	0.350
	0.242
	0.338
	0.463
	0.352
	0.320
	0.333
	0.399
	0.291
	0.364
	0.312
	0.233
	0.298
	0.165
	0.300
	0.380
	0.318
	0.257
	0.230
	0.185
	0.155
	0.615
	0.667
	0.675
	0.664
	0.598
	0.607
	0.597
	0.473
	0.483
	0.583
	0.564
	0.754
	0.779
	0.740
	0.694
];

%}



%% ========================================================================
% Analyze the NC 13 nuclear concentration, for normalization purposes
% =========================================================================
% The reason why we're normalizing to nc13 is because all time courses have
% a complete nc13 and there's not that much variation within nc13, so
% taking the average should be sufficient (ie, we probably don't have to do
% a curve-fitting procedure).
% {

%
% Calc the nc 13 average
%
n_embryos = length(Soln);
Cnuc13 = NaN(1,n_embryos);
NCRI13 = NaN(1,n_embryos);
for iii = 1:n_embryos
	t1 = Soln(iii).t;
	Cnuc1 = 1./(V_eff*Soln(iii).Anuc)/ 0.602;
	Inuc1 = Soln(iii).nucsignal;
	Icyt1 = Soln(iii).cytsignal;

	%
	% Get anchorpoint, which is the beginning of nc14
	%
	t14_0 = t1(1,end);

	%
	% Get time points for nc13
	%
	t_1 = t1(:,4); 
	if all(isnat(t_1))
		continue
	end
	tplot = minutes(t_1-t14_0);
	[tplot,iplot] = sort(tplot);
	vnan = isnan(tplot);
	tplot(vnan) = [];
	iplot(vnan) = [];
	m = (tendmean(4) - tendmean(3))/(tend(iii,4) - tend(iii,3));
	t_hat = m*(tplot - tend(iii,4)) + tendmean(4);

	%
	% Get nuc Bcd conc and NCRI for nc13, and average over only the final
	% 5 minutes
	%
	Cnuc_1 = Cnuc1(iplot,4);
	Inuc_1 = Inuc1(iplot,4);
	Icyt_1 = Icyt1(iplot,4);
	v = t_hat(end)-t_hat < 5;
	Cnuc13(iii) = mean(Cnuc_1(v),'omitmissing');
	NCRI13(iii) = mean(Inuc_1(v)./Icyt_1(v),'omitmissing');

end
Cnuc13bar = mean(Cnuc13,'omitmissing');
NCRI13bar = mean(NCRI13,'omitmissing');


%
% NCRI13 is very tightly linear wrt APcoord. We will use this relationship
% to fill in the NaN's of the APcoord
%
v = ~isnan(s) & ~isnan(NCRI13');
[m,b,stats] = linlsq(s(v),NCRI13(v)');
v = isnan(s);
s(v) = (NCRI13(v) - b)/m;


%}


%% ========================================================================
% Create a series of time course plots for each bin in "s"
% =========================================================================
% {
% The full timecourse embryos are the first 25 elements of Soln, and are
% plotted only up to AP coordinate of 0.50
Soln0 = Soln(1:25);
s0 = s(1:25);


%
% Create bins and mesh
%
s1 = 0.13;
s2 = 0.48;
ds = 0.07;
s_mesh = (s1:ds:s2)';
sb1 = max(0,s_mesh - ds/2);
sb2 = s_mesh + ds/2; sb2(end) = 1;
ns = length(sb1);

sem_flag = false; 
yesplot = false;

%
% Bin time course embryos into the stated mesh from AP = 0.15 to 0.50
%
n = zeros(ns,1);
for i = 1:ns
	v = s0 >= sb1(i) & s0 < sb2(i);
	if sum(v) <= 1
		v = false(size(v));
	end
	n(i) = sum(v);
	
	data_bcd_timecourse(i) = extract_RICS_timecourse(Soln0(v),yesplot);
end

%
% Make plots
%
figure('pos',[9         102        1899         592])
TL = tiledlayout(2,5,'TileSpacing','Compact','Padding','Compact');
C = colormap('lines');
t_end = 40; % shift zero to the final time point (coincides with gastrulation)

for i = 1:ns

	%
	% Extract variables
	%
	struct2vars(data_bcd_timecourse(i))
	t = t - t_end;

	%
	% First row
	%
	nexttile(1)
	errorbar(t(:),Cnuc(:),Cnuc_e(:))
	hold on
	xlabel('time [min]')
	ylabel('tot conc [nM]')
	YLIM = ylim;
	ylim([0 YLIM(2)])

	nexttile(2)
	errorbar(t(:),Cfree(:),Cfree_e(:))
	hold on
	xlabel('time [min]')
	ylabel('free conc [nM]')
	YLIM = ylim;
	ylim([0 YLIM(2)])

	nexttile(3)
	errorbar(t(:),CDNAbound(:),CDNAbound_e(:))
	hold on
	xlabel('time [min]')
	ylabel('correl conc [nM]')
	YLIM = ylim;
	ylim([0 YLIM(2)])

	nexttile(4)
	errorbar(t(:),Cuncorrel(:),Cuncorrel_e(:))
	hold on
	xlabel('time [min]')
	ylabel('uncorrel conc [nM]')
	YLIM = ylim;
	ylim([0 YLIM(2)])
	
	% For the dose/response plots, ignore the first handful of points in each
	% nc (pre-pseudo steady state). This is achieved by the array "idx_de"
	nexttile(5)
	errorbar(Cfree(idx_de),CDNAbound(idx_de),CDNAbound_e(idx_de),CDNAbound_e(idx_de),Cfree_e(idx_de),Cfree_e(idx_de))
	hold on
	xlabel('free conc [nM]')
	ylabel('correl conc [nM]')
	YLIM = ylim;
	ylim([0 YLIM(2)])

	%
	% second row
	%
	nexttile(6)
	errorbar(t(:),Cnuc_ch(:),Cnuc_ch_e(:))
	hold on
	xlabel('time [min]')
	ylabel('H2Av conc [nM]')
	YLIM = ylim;
	ylim([0 YLIM(2)])

	nexttile(7)
	errorbar(t(:),phicc(:),phicc_e(:))
	hold on
	xlabel('time [min]')
	ylabel('free frac [nM]')
	ylim([0 1])

	nexttile(8)
	errorbar(t(:),phicc(:),phicc_e(:))
	hold on
	xlabel('time [min]')
	ylabel('correl frac [nM]')
	ylim([0 1])

	nexttile(9)
	errorbar(t(:),phiuncorrel(:),phiuncorrel_e(:))
	hold on
	xlabel('time [min]')
	ylabel('uncorrel frac [nM]')
	ylim([0 1])

	% dose/response
	nexttile(10)
	errorbar(Cfree(idx_de),Cuncorrel(idx_de),Cuncorrel_e(idx_de),Cuncorrel_e(idx_de),Cfree_e(idx_de),Cfree_e(idx_de))
	hold on
	xlabel('free conc [nM]')
	ylabel('uncorrel conc [nM]')
	YLIM = ylim;
	ylim([0 YLIM(2)])


end



save Mat/plotoutput_Bcd_orientation data_bcd_timecourse
%}




%% ========================================================================
% Run the plotting function for nc14 at various AP coordinates
% =========================================================================
% {

yesplot = true;
data_bcd_orientation = extract_RICS_timecourse_orientation(Soln,s,yesplot);

save Mat/plotoutput_Bcd_orientation data_bcd_orientation
%}















