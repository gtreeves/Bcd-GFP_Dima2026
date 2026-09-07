function data = plot_RICS_timecourse(Soln,outfilename,normalize,convolve)

if ~exist('outfilename','var') || ~(ischar(outfilename) || isstring(outfilename))
    outfilename = '';
else
    outfilename = char(outfilename);
end

if ~isempty(outfilename) && ~strcmp(outfilename(end),'_')
	outfilename(end+1) = '_';
end

if ~exist('normalize','var') || (~isnumeric(normalize) && ~islogical(normalize))
	normalize = false;
elseif isnumeric(normalize)
	warning('"normalize" should be a logical variable. Converting to logical')
	normalize = ~~normalize;
end

if ~exist('convolve','var') || (~isnumeric(convolve) && ~islogical(convolve))
	convolve = false;
elseif isnumeric(convolve)
	warning('"convolve" should be a logical variable. Converting to logical')
	convolve = ~~convolve;
end

% %
% % Histone controls
% %
% load('C:\Users\gtreeves\Dropbox\Matlab\RICS\control_RFP\Mat\average_histone','t','Cnucbar','Cnucstd')
% t_H2A = t;
% H2Abar = Cnucbar;
% H2Astd = Cnucstd;

%
% Load timecourses
%
if ischar(Soln) && exist(Soln,'dir')
	load(Soln,'Soln')
end

%
% Geometric params
%
w00 = Soln(1).metadata{end}{1}.w0;
w0r = 0.3;
chi = 3;
gamma = sqrt(2)/4; % Brown2008 has this as sqrt(2)/4, not sqrt(2)/2
% I am not sure where I got the "sqrt(2)/2" from. So I changed it.
V_eff = pi^1.5*chi*w00^3/gamma;
Vnuc_ch = pi^1.5*chi*w0r^3/gamma;
w0_cc = sqrt(0.5*(w00.^2 + w0r.^2));
Vcc = pi^1.5*chi*w0_cc.^3/gamma;


%% ========================================================================
% Normalize nuclear cycle lengths
%
% So, it turns out that these embryos all had differences in their nuclear
% cycle length. So we will have to analyze each of them and normalize their
% lengths to match each other, similar to what was done in Reeves et al.,
% 2012. 
%
% =========================================================================
% {


NC_list = 10:14;
n_embryos = length(Soln);
tstart = NaN(n_embryos,length(NC_list)); % array for start times of different NCs for different embryos
tend = NaN(n_embryos,length(NC_list)); % array for end times of different NCs for different embryos
for iii = 1:n_embryos
	t1 = Soln(iii).t;
	t14_0 = min(t1(:,end)); % start of nc 14

	tstart(iii,:) = minutes(min(t1) - t14_0);
	tend(iii,:) = minutes(max(t1) - t14_0);

end

%
% Create general mesh of time points, fixed
%
dt = 2; % roughly two minutes between each time point
dti = [2 4 6 12 40]; % used to be 40. % the desired durations for nc10 through nc14
dtm = [6 6 6 8];%the desired gaps between the starts of adjacent cycles
tstartmean = [0 cumsum(dti(1:end-1)+dtm)];
tendmean = tstartmean + dti;
t14mesh0 = tstartmean(end);
tstartmean = tstartmean - t14mesh0;
tendmean = tendmean - t14mesh0;
Delta_t = dti;
nt = round(Delta_t/dt) + 1;
t_mesh = NaN(max(nt),length(NC_list));
for i = 1:length(NC_list)
	t_mesh(1:nt(i),i) = linspace(tstartmean(i),tendmean(i),nt(i))';
end

% t_mesh

%}


%% ========================================================================
% Do a convolution of the nuc_ch in nc14
% =========================================================================
% {

if convolve

	%
	% To do the convolution, we first have to extract the nuc_ch in nc14 for
	% each embryo and put them onto the same mesh.
	%
	Cbar = zeros(max(nt),n_embryos);
	for iii = 1:n_embryos
		t1 = Soln(iii).t;
		t_1 = t1(:,end); t_1 = t_1(:);
		t14_0 = min(t1(:,end)); % start of nc 14
		% Cnuc1 = 1./(V_eff*Soln(iii).Anuc)/ 0.602;
		Cnuc_ch1 = 1./(Vnuc_ch*Soln(iii).Anuc_ch)/ 0.602;
		Cnuc_ch1(Soln(iii).Anuc_ch < 1e-5) = NaN;

		

		%
		% Take the time points in the current NC and normalize them wrt
		% to the end of nc14 based on the average nc duration.
		%
		tplot = minutes(t_1-t14_0);
		[tplot,iplot] = sort(tplot);
		vnan = isnan(tplot);
		tplot(vnan) = [];
		iplot(vnan) = [];

		m = (tendmean(4) - tendmean(3))/(tend(iii,4) - tend(iii,3));
		t_hat = m*(tplot - tstart(iii,end)) + tstartmean(end);

		%
		% Load variables
		%
		C1 = Cnuc_ch1(:,end);
		C1 = C1(:);
		C1 = C1(iplot);

		%
		% Building averages for some of the variables
		%
		C1 = interp1(t_hat,C1,t_mesh(:,end));
		Cbar(:,iii) = C1;

	end

	t1 = t_mesh(:,end);

	%
	% Testing to see if a convolution in the nuc_ch will help us figure out
	% what time point each embryo started their nc 14 curve with. It looks
	% like embryo 2 starts the earliest, so all others will be compared to
	% that one, even though it also ends too early.

    %AGP: This earlier comment was from Greg using Sadia's Zld/GAF RICS
    %data, I think I will have to change it for CRISPR-Med-GFP data.

    %Update:Not using convolve for our data. But in the future if I am need
    %to think about this.
	%
	J2 = Cbar(:,2);
	v = isnan(J2);
	J2(v) = 0;
	J2 = fft(J2);

	G1 = zeros(length(t1),n_embryos);
	for iii = 1:n_embryos
		J1 = Cbar(:,iii);
		v = isnan(J1);
		J1(v) = 0;
		J1 = fft(J1);
		G1(:,iii) = ifft(J2.*conj(J1));
	end
	[~,imax] = max(G1);

else
	imax = ones(1,n_embryos);

end

%}



%% ========================================================================
% Make averages of data for all embryos, with some adjustments
% =========================================================================
% {

%
% Interpolate arrays onto the standard t_mesh
%
Cnuc_interp = NaN(max(nt),length(NC_list),n_embryos);
Ccyt_interp = NaN(max(nt),length(NC_list),n_embryos);
Cnuc_ch_interp = NaN(max(nt),length(NC_list),n_embryos);
phinuc_interp = NaN(max(nt),length(NC_list),n_embryos);
phicyt_interp = NaN(max(nt),length(NC_list),n_embryos);
phicc_interp = NaN(max(nt),length(NC_list),n_embryos);
Inuc_interp = NaN(max(nt),length(NC_list),n_embryos);
Icyt_interp = NaN(max(nt),length(NC_list),n_embryos);
for iii = 1:n_embryos

	t1 = Soln(iii).t;
	t14_0 = min(t1(:,end));

	Cnuc1 = 1./(V_eff*Soln(iii).Anuc)/ 0.602;
	Cnuc1(Soln(iii).Anuc < 1e-5) = NaN;
	Ccyt1 = 1./(V_eff*Soln(iii).Acyt)/ 0.602;
	Ccyt1(Soln(iii).Acyt < 1e-5) = NaN;
	Cnuc_ch1 = 1./(Vnuc_ch*Soln(iii).Anuc_ch)/ 0.602;
	Cnuc_ch1(Soln(iii).Anuc_ch < 1e-5) = NaN;
	phinuc1 = Soln(iii).phinuc;
	phicyt1 = Soln(iii).phicyt;
	phicc1 = Vcc./Vnuc_ch.*Soln(iii).Acc./Soln(iii).Anuc_ch;
    % Agp: Adding this if statement for our purposes, might not be needed
    % in other cases, cause correlated concentration cannot be less than
    % zero
    phicc1(Soln(iii).Acc < 0) = 0;
	Inuc1 = Soln(iii).nucsignal;
	Icyt1 = Soln(iii).cytsignal;

	%
	% Process the stages
	%
	% nfiles = length(Cnuc1);
	for i = 1:length(NC_list)
		
		%
		% Take the time points in the current NC and normalize them wrt
		% to the end of nc14 based on the average nc duration.
		%
		t_1 = t1(:,i); 
		if all(isnat(t_1)) % break out of NC's that have no data
			continue
		end
		tplot = minutes(t_1-t14_0); %t_plot is the original embryo-specific time axis
		[tplot,iplot] = sort(tplot);
		vnan = isnan(tplot);
		tplot(vnan) = [];
		iplot(vnan) = [];

		if i < length(NC_list)
			m = (tendmean(i) - tstartmean(i))/(tend(iii,i) - tstart(iii,i));

		% special cases for nc 14, because not all gastrulation captured:
		elseif ~isnan(tend(iii,4)) && ~isnan(tend(iii,3))
			% This case assumes there exist data in both nc 12 and nc
			% 13...this could be a problem because some embryos have only
			% nc 13,14, and others only nc14. In that case, we fall back on
			% alternatives.
			m = (tendmean(4) - tendmean(3))/(tend(iii,4) - tend(iii,3)); % scaling factor of nc 13 to nc 12 applied to NC 14

		elseif ~isnan(tend(iii,4))
			% The case for embryos that have nc13 but not 12:
			m = (tstartmean(5) - tstartmean(4))/(tstart(iii,5) - tstart(iii,4));

		else
			% final alternative: if there's only nc14...don't change
			% anything
			m = 1;			
		end
		t_hat = m*(tplot - tstart(iii,i)) + tstartmean(i); %t_hat is the normalized time axis

		%
		% Time shift for nc 14
		%
		if i == length(NC_list)
			t0 = t_mesh(imax(iii),end);
		else
			t0 = 0;
		end
		if all(isnan(t_hat)) % break out of NC's that have no data
			continue
		end

		%
		% Do the interpolation. If only one data point, no interpolation
		% needed
		%
		if length(t_hat) > 1
			Cnuc_interp(1:nt(i),i,iii) = interp1(t_hat+t0,Cnuc1(iplot,i),t_mesh(1:nt(i),i));
			Ccyt_interp(1:nt(i),i,iii) = interp1(t_hat+t0,Ccyt1(iplot,i),t_mesh(1:nt(i),i));
			Cnuc_ch_interp(1:nt(i),i,iii) = interp1(t_hat+t0,Cnuc_ch1(iplot,i),t_mesh(1:nt(i),i));
			phinuc_interp(1:nt(i),i,iii) = interp1(t_hat+t0,phinuc1(iplot,i),t_mesh(1:nt(i),i));
			phicyt_interp(1:nt(i),i,iii) = interp1(t_hat+t0,phicyt1(iplot,i),t_mesh(1:nt(i),i));
			phicc_interp(1:nt(i),i,iii) = interp1(t_hat+t0,phicc1(iplot,i),t_mesh(1:nt(i),i));
			Inuc_interp(1:nt(i),i,iii) = interp1(t_hat+t0,Inuc1(iplot,i),t_mesh(1:nt(i),i));
			Icyt_interp(1:nt(i),i,iii) = interp1(t_hat+t0,Icyt1(iplot,i),t_mesh(1:nt(i),i));
		else
			Cnuc_interp(1:nt(i),i,iii) = Cnuc1(iplot,i);
			Ccyt_interp(1:nt(i),i,iii) = Ccyt1(iplot,i);
			Cnuc_ch_interp(1:nt(i),i,iii) = Cnuc_ch1(iplot,i);
			phinuc_interp(1:nt(i),i,iii) = phinuc1(iplot,i);
			phicyt_interp(1:nt(i),i,iii) = phicyt1(iplot,i);
			phicc_interp(1:nt(i),i,iii) = phicc1(iplot,i);
			Inuc_interp(1:nt(i),i,iii) = Inuc1(iplot,i);
			Icyt_interp(1:nt(i),i,iii) = Icyt1(iplot,i);

		end
	end
end

% %
% % Correct phicc just within nc 14
% %
% if i == length(NC_list)
% 	Cnuc_ch0 = interp1(t_H2A2{5},H2A_out{5}/2,t_hat+t0);
% 	phiccplot = phiccplot.*Cnuc_chplot./Cnuc_ch0;
% end

NCRI = Inuc_interp./Icyt_interp;

%
% % Delete embryos with NCRI < 1 and the one aberrant embryo (for Dorsal),
% % and also ones that are full NaN
% % AGP: Not deleting embryos with NCRI<1 anymore
% NCRI1 = meanDU(squeeze(meanDU(NCRI)))';
% Cnuc13 = meanDU(squeeze(Cnuc_interp(:,4,:)))';
NCRI1  = meanDU(meanDU(NCRI,1),2);          % 1 x 1 x n_embryos
NCRI1  = reshape(NCRI1,1,[]);               % 1 x n_embryos

Cnuc13 = meanDU(Cnuc_interp(:,4,:),1);      % 1 x 1 x n_embryos
Cnuc13 = reshape(Cnuc13,1,[]);              % 1 x n_embryos
% if sum(Cnuc13 > 140) == 1 && n_embryos > 5
% 	v = NCRI1 < 1 | Cnuc13 > 140;
% else
% 	v = NCRI1 < 1;
% end


% v_allnan = all(all(isnan(Cnuc_interp),1),2); % 1 x 1 x n_embryos
% v_allnan = reshape(v_allnan,1,[]);           % 1 x n_embryos

v = all(squeeze(all(isnan(Cnuc_interp))))';
Cnuc_interp(:,:,v) = [];
Ccyt_interp(:,:,v) = [];
Cnuc_ch_interp(:,:,v) = [];
phinuc_interp(:,:,v) = [];
phicyt_interp(:,:,v) = [];
phicc_interp(:,:,v) = [];
Inuc_interp(:,:,v) = [];
Icyt_interp(:,:,v) = [];
NCRI(:,:,v) = [];
Cnuc13(v) = [];

n_embryos = size(Cnuc_interp,3);


t = [t_mesh;NaN(1,length(NC_list))]; % add extra NaN to separate NCs, for plotting purposes,
                                    % when we do reshape and stack the columns NaN is breakpoint
t0 = -max(t(:,end));


%Agp : What is this? How to change?
i1 = [3 3 4 5 8];
% i1 = [2 2 3 4 5];
idx_de = true(size(t));
% for i = 1:length(NC_list)
% 	idx_de(i1(i):end,i) = true;
% end


if normalize
	chi0 = meanDU(Cnuc13)./Cnuc13; % need to think about dimensions before doing this
                                   % Agp: What about embryos which don't have NC-13 data do we not normalize them  

	chi = repmat(permute(chi0,[2,3,1]),size(t,1),length(NC_list),1);
	chi = reshape(chi,[(max(nt)+1)*length(NC_list),n_embryos]);
else
	chi = 1;
end


%agp: Debugging thing only for understanding interpolation
%tplot_at_mesh = interp1(t_hat, t_plot, t_mesh(7,4));


%
% Create variables for plotting
%
Cnuc = [Cnuc_interp;NaN(1,length(NC_list),n_embryos)]; % pad bottom of each nc with a NaN
Cnucplot = reshape(Cnuc,[(max(nt)+1)*length(NC_list),n_embryos]);
Cnucbar = meanDU(Cnuc,3);
Cnucstd = semDU(Cnuc,3);

Ccyt = [Ccyt_interp;NaN(1,length(NC_list),n_embryos)];
Ccytplot = reshape(Ccyt,[(max(nt)+1)*length(NC_list),n_embryos]);
Ccytbar = meanDU(Ccyt,3);
Ccytstd = semDU(Ccyt,2);

Cnuc_ch = [Cnuc_ch_interp;NaN(1,length(NC_list),n_embryos)];
Cnuc_chplot = reshape(Cnuc_ch,[(max(nt)+1)*length(NC_list),n_embryos]);
Cnuc_chbar = meanDU(Cnuc_ch,3);
Cnuc_chstd = semDU(Cnuc_ch,3);

phinuc = [phinuc_interp;NaN(1,length(NC_list),n_embryos)];
phinucplot = reshape(phinuc,[(max(nt)+1)*length(NC_list),n_embryos]);
phinucbar = meanDU(phinuc,3);
phinucstd = semDU(phinuc,3);

phicyt = [phicyt_interp;NaN(1,length(NC_list),n_embryos)];
phicytplot = reshape(phicyt,[(max(nt)+1)*length(NC_list),n_embryos]);
phicytbar = meanDU(phicyt,3);
phicytstd = semDU(phicyt,3);

phicc = [phicc_interp;NaN(1,length(NC_list),n_embryos)];
phiccplot = reshape(phicc,[(max(nt)+1)*length(NC_list),n_embryos]);
phiccbar = meanDU(phicc,3);
phiccstd = semDU(phicc,3);

NCRI = [NCRI;NaN(1,length(NC_list),n_embryos)];
NCRIplot = reshape(NCRI,[(max(nt)+1)*length(NC_list),n_embryos]);
NCRIbar = meanDU(NCRI,3);
NCRIstd = semDU(NCRI,3);

Inuc = [Inuc_interp;NaN(1,length(NC_list),n_embryos)]; % pad bottom of each nc with a NaN
Inucplot = reshape(Inuc,[(max(nt)+1)*length(NC_list),n_embryos]);
Inucbar = meanDU(Inuc,3);
Inucstd = semDU(Inuc,3);

Icyt = [Icyt_interp;NaN(1,length(NC_list),n_embryos)];
Icytplot = reshape(Icyt,[(max(nt)+1)*length(NC_list),n_embryos]);
Icytbar = meanDU(Icyt,3);
Icytstd = semDU(Icyt,3);


Cfree = Cnuc.*(1-phinuc);
Cfreeplot = Cnucplot.*(1-phinucplot);
Cfreebar = Cnucbar.*(1-phinucbar);
Cfreestd = Cfreebar.*sqrt((Cnucstd./Cnucbar).^2 + (phinucstd./phinucbar).^2);

CDNAbound = Cnuc.*phicc;
CDNAboundplot = Cnucplot.*phiccplot;
CDNAboundbar = Cnucbar.*phiccbar;
CDNAboundstd = CDNAboundbar.*sqrt((Cnucstd./Cnucbar).^2 + (phiccstd./phiccbar).^2);

Cuncorrel = Cnuc.*(phinuc - phicc);
Cuncorrelplot = Cnucplot.*(phinucplot - phiccplot);
Cuncorrelbar = Cnucbar.*(phinucbar - phiccbar);
Cuncorrelstd = Cuncorrelbar.*sqrt((Cnucstd./Cnucbar).^2 + (phinucstd.^2 + phiccstd.^2)./(phinucbar-phiccbar).^2);


%}

%% ========================================================================
% Make plots of data for all embryos, with some adjustments
% =========================================================================
% {

%
% Make a plot of 10 tiles. 
%
figure('pos',[9         102        1899         592])
TL = tiledlayout(2,5,'TileSpacing','Compact','Padding','Compact');


%
% Make cell variables for each of the aplots (X,Y,S)
%
NAN = NaN(size(Cfreebar));
NAN(idx_de) = 1;

X = [repmat({t(:)+t0},1,4),Cfreeplot.*chi
	repmat({t(:)+t0},1,4),Cfreeplot.*chi]; % the "1" is a placeholder since the final plot is different from the others
Y = {Cnucplot.*chi, Cfreeplot.*chi, CDNAboundplot.*chi, Cuncorrelplot.*chi, CDNAboundplot.*chi
	Cnuc_chplot.*chi, 1-phinucplot, phiccplot, phinucplot-phiccplot, Cuncorrelplot.*chi};


Xe = [repmat({t(:)+t0},1,4),Cfreebar(idx_de(:))
	repmat({t(:)+t0},1,4),Cfreebar(idx_de(:))]; % the "1" is a placeholder since the final plot is different from the others
Ye = {Cnucbar(:), Cfreebar(:), CDNAboundbar(:), Cuncorrelbar(:), CDNAboundbar(idx_de(:))
	Cnuc_chbar(:), 1-phinucbar(:), phiccbar(:), phinucbar(:)-phiccbar(:), Cuncorrelbar(idx_de(:))};
Sx = {0, 0, 0, 0, Cfreestd(idx_de(:))
	0, 0, 0, 0, Cfreestd(idx_de(:))};
Sy = {Cnucstd(:), Cfreestd(:), CDNAboundstd(:), Cuncorrelstd(:), CDNAboundstd(idx_de(:))
	Cnuc_chstd(:), phinucstd(:), phiccstd(:), sqrt(phinucstd(:).^2+phiccstd(:).^2), Cuncorrelstd(idx_de(:))};

%
% Make graph-labeling variables
%
Xlabel = [repmat({'Time until gastrulation [min]'},1,4),{'Free nuc conc'}
	repmat({'Time until gastrulation [min]'},1,4),{'Free nuc conc'}];
Ylabel = {'Total nuc conc. [nM]', 'Free nuc conc. [nM]', 'Correlated conc. [nM]', 'Uncorrelated conc. [nM]', 'Correlated conc.'
	'Histone conc [nM]', 'free (1-\phi_{2c,nuc})', '\phi_{cc}', 'uncorrel (\phi_{2c,nuc}-\phi_{cc})', 'Uncorrelated conc'};
YLIM2 = [inf(1,5)
	inf,1,1,1,inf];


% Double for loop: Plot two rows of 5 tiles each

for i = 1:2
	for j = 1:5
	nexttile % 1,1
	plot(X{i,j},Y{i,j},'o-','linewidth',0.5)
	hold on
    % Adding error bars 
	if j < 5
		errorbar(Xe{i,j},Ye{i,j},Sy{i,j},'o-k','linewidth',2)
	else
		errorbar(Xe{i,j},Ye{i,j},Sy{i,j},Sy{i,j},Sx{i,j},Sx{i,j},'o-k','linewidth',2)
		xlim([0 inf])
	end
	xlabel(Xlabel{i,j})
	ylabel(Ylabel{i,j})
	ylim([0 YLIM2(i,j)])
	set(gca,'fontsize',12)

	end
end
lgd = legend("Embryo " + string(1:n_embryos));
lgd.Layout.Tile = 'east';
% disp(n_embryos)
% disp(size(Cnucplot))
% disp(outfilename)
disp(['Figs/',outfilename,'allembryos'])
epsDU(gcf,['Figs/',outfilename,'allembryos'],12)
% save Mat/average_Zld t Cnucbar Cnucstd


%}



%% ========================================================================
% Make plot of average curves
% =========================================================================
% {

%
% Make a plot of 10 tiles. 
%
figure('pos',[9         102        1899         592])
    TL = tiledlayout(2,5,'TileSpacing','Compact','Padding','Compact');

%
% Make cell variables for each of the plots (X,Y,S)
%
NAN = NaN(size(Cfreebar));
NAN(idx_de) = 1;

X = [repmat({t+t0},1,4),Cfreebar.*NAN
	repmat({t+t0},1,4),Cfreebar.*NAN]; % the "1" is a placeholder since the final plot is different from the others
Y = {Cnucbar, Cfreebar, CDNAboundbar, Cuncorrelbar, CDNAboundbar.*NAN
	Cnuc_chbar, 1-phinucbar, phiccbar, phinucbar-phiccbar, Cuncorrelbar.*NAN};
Sx = {0, 0, 0, 0, Cfreestd.*NAN
	0, 0, 0, 0, Cfreestd.*NAN};
Sy = {Cnucstd, Cfreestd, CDNAboundstd, Cuncorrelstd, CDNAboundstd.*NAN
	Cnuc_chstd, phinucstd, phiccstd, sqrt(phinucstd.^2+phiccstd.^2), Cuncorrelstd.*NAN};

%
% Make graph-labeling variables
%
Xlabel = [repmat({'time until gast. [min]'},1,4),{'Free nuc conc'}
	repmat({'time until gast. [min]'},1,4),{'Free nuc conc'}];
Ylabel = {'Total nuc conc. [nM]', 'Free nuc conc. [nM]', 'Correlated conc. [nM]', 'Uncorrelated conc. [nM]', 'Correlated conc'
	'Histone conc [nM]', 'free (1-\phi_{2c,nuc})', '\phi_{cc}', 'uncorrel (\phi_{2c,nuc}-\phi_{cc})', 'Uncorrelated conc'};
YLIM2 = [inf(1,5)
inf,1,1,1,inf];

%
% Double for loop: Plot two rows of 5 tiles each
%
for i = 1:2
	for j = 1:5
	nexttile % 1,1

	if j < 5
		errorbar(X{i,j},Y{i,j},Sy{i,j},'o-','linewidth',2)
	else
		errorbar(X{i,j},Y{i,j},Sy{i,j},Sy{i,j},Sx{i,j},Sx{i,j},'o-','linewidth',2)
		xlim([0 inf])
	end
	xlabel(Xlabel{i,j})
	ylabel(Ylabel{i,j})
	ylim([0 YLIM2(i,j)])
	set(gca,'fontsize',12)

	end
end
lgd = legend("NC " + string(NC_list));
lgd.Layout.Tile = 'east';

epsDU(gcf,['Figs/',outfilename,'avgembryo'],12)

%}



%% ========================================================================
% Make plot of intensities
% =========================================================================
% {

%
% Make a plot of 10 tiles. 
%
figure('pos',[9         102        1899         592])
TL = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
C = colormap('lines');

%
% Make cell variables for each of the plots (X,Y,S)
%
NAN = NaN(size(Cfreebar));
NAN(idx_de) = 1;

X = [repmat({t+t0},1,3)];
Y = {Inucbar, Icytbar, NCRIbar};
Sx = {0, 0, 0};
Sy = {Inucstd, Icytstd, NCRIstd};

%
% Make graph-labeling variables
%
Xlabel = [repmat({'time until gast. [min]'},1,3)];
Ylabel = {'Nuc intensity [AU]', 'Cyt intensity [AU]','NCRI'};
YLIM2 = [inf(1,3)];

%
% Double for loop: Plot two rows of 5 tiles each
%
for i = 1:1
	for j = 1:3
	nexttile % 1,1

	if j < 5
		errorbar(X{i,j},Y{i,j},Sy{i,j},'o-','linewidth',2)
	else
		errorbar(X{i,j},Y{i,j},Sy{i,j},Sy{i,j},Sx{i,j},Sx{i,j},'o-','linewidth',2)
		xlim([0 inf])
	end
	xlabel(Xlabel{i,j})
	ylabel(Ylabel{i,j})
	ylim([0 YLIM2(i,j)])
	set(gca,'fontsize',12)

	end
end

lgd = legend("NC " + string(NC_list));
lgd.Layout.Tile = 'east';
epsDU(gcf,['Figs/',outfilename,'intensities'],12)

%}





%% ========================================================================
% Make plots of correlated and uncorrelated conc vs. Total nuc conc
% =========================================================================
% {
figure
x = Cnucbar; y = CDNAboundbar; xs = Cnucstd; ys = CDNAboundstd;
x(~idx_de) = NaN; y(~idx_de) = NaN; xs(~idx_de) = NaN; ys(~idx_de) = NaN;
errorbar(x,y,ys,ys,xs,xs,'o-','linewidth',2)
xlabel('Total nuc conc')
ylabel('Correlated conc')
xlim([0 inf])
ylim([0 inf])
set(gca,'fontsize',12)
lgd = legend("NC " + string(NC_list));
% lgd.Layout.Tile = 'east';


figure
for i = 1:length(NC_list)
	h = errorbar(Cnucbar(i1(i):end,i),Cuncorrelbar(i1(i):end,i),...
		Cuncorrelstd(i1(i):end,i),Cuncorrelstd(i1(i):end,i),Cnucstd(i1(i):end,i),Cnucstd(i1(i):end,i),'o-');
	hold on
	for j = 1:length(h)
		set(h(j),'linewidth',2)
	end
end
% errorbar(Cfreebar(idx_de(:)),Cuncorrelbar(idx_de(:)),...
% 	Cuncorrelstd(idx_de(:)),Cuncorrelstd(idx_de(:)),...
% 	Cfreestd(idx_de(:)),Cfreestd(idx_de(:)),'o-','linewidth',2)
xlim([0 inf])
ylim([0 inf])
xlabel('Total nuc conc')
ylabel('Uncorrelated conc')
set(gca,'fontsize',12)

epsDU(gcf,['Figs/',outfilename,'totalnucconc_vs_uncorrelconc'],12)

%}

%% ========================================================================
% Save output variables
% =========================================================================
% {

data.t = t;
data.t0 = t0;
data.idx_de = idx_de;
data.Cnuc = Cnucbar;
data.Ccyt = Ccytbar;
data.Cfree = Cfreebar;
data.CDNAbound = CDNAboundbar;
data.Cuncorrel = Cuncorrelbar;
data.Cnuc_ch = Cnuc_chbar;
data.phinuc = phinucbar;
data.phicyt = phicytbar;
data.phicc = phiccbar;
data.phiuncorrel = phinucbar-phiccbar;
data.NCRI = NCRIbar;
data.Inuc = Inucbar;
data.Icyt = Icytbar;

data.Cnuc_e = Cnucstd;
data.Ccyt_e = Ccytstd;
data.Cfree_e = Cfreestd;
data.CDNAbound_e = CDNAboundstd;
data.Cuncorrel_e = Cuncorrelstd;
data.Cnuc_ch_e = Cnuc_chstd;
data.phinuc_e = phinucstd;
data.phicyt_e = phicytstd;
data.phicc_e = phiccstd;
data.phiuncorrel_e = sqrt(phinucstd.^2+phiccstd.^2);
data.NCRI_e = NCRIstd;
data.Inuc_e = Inucstd;
data.Icyt_e = Icytstd;

data.full.Cnuc = Cnuc;
data.full.Ccyt = Ccyt;
data.full.Cfree = Cfree;
data.full.CDNAbound = CDNAbound;
data.full.Cnuc_ch = Cnuc_ch;
data.full.phinuc = phinuc;
data.full.phicyt = phicyt;
data.full.phicc = phicc;
data.full.NCRI = NCRI;
data.full.Inuc = Inuc;
data.full.Icyt = Icyt;
data.full.chi = chi;



%}




