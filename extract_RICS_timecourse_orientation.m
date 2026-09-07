function data = extract_RICS_timecourse_orientation(Soln,s,yesplot,sem_flag)


%
% Load timecourses if "Soln" is passed as a path to a properly-formatted
% mat file containing a structure "Soln", which is the output of our RICS
% timecourse analysis pipeline (but the usual expectation is that "Soln" is
% the structure itself).
%
if ischar(Soln) && exist(Soln,'dir')
	load(Soln,'Soln')
end

%
% Flag to determine whether we use the SEM or STD for the errobar
%
if ~exist('sem_flag','var') || isempty(sem_flag) || ~isscalar(sem_flag)
	sem_flag = true; % therefore, use SEM (default)
end

%
% Geometric params
%
w00 = Soln(1).metadata{end}{1}.w0;
w0r = 0.3;
chi = 3;
gamma = sqrt(2)/4; 
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
n_ncs = length(NC_list);
n_embryos = length(Soln);
tstart = NaN(n_embryos,n_ncs);
tend = NaN(n_embryos,n_ncs);
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
dti = [2 4 6 12 40]; % used to be 40.
dtm = [6 6 6 8];
tstartmean = [0 cumsum(dti(1:end-1)+dtm)];
tendmean = tstartmean + dti;
t14mesh0 = tstartmean(end);
tstartmean = tstartmean - t14mesh0;
tendmean = tendmean - t14mesh0;
Delta_t = dti;
nt = round(Delta_t/dt) + 1;
t_mesh = NaN(max(nt),n_ncs);
for i = 1:n_ncs
	t_mesh(1:nt(i),i) = linspace(tstartmean(i),tendmean(i),nt(i))';
end

t = [t_mesh;NaN(1,n_ncs)]; % add extra NaN to separate NCs
t0 = -max(t(:,end));


% Create "start" vector indicating the indices for each nuclear cycle
% before which oddities of exiting mitosis make the data empirically
% "weird"
i1 = [3 3 4 5 8];
idx_de = false(size(t));
for i = 1:n_ncs
	idx_de(i1(i):end,i) = true;
end

%}



%% ========================================================================
% Interpolate basic data for all embryos onto the same time mesh
% =========================================================================
% Interpolate basic variables onto the standard t_mesh. Note the length of
% the first dimension (time) is max(nt)+1. The "plus one" is to pad the
% bottom with a NaN.
%
% By "basic" I mean Inuc, Icyt, Varnuc, Varcyt, Anuc, Acyt, Anuc_ch, Acc,
% phinuc, phicyt, etc. Full list just below:
% {

varnames = {'Inuc', 'Icyt', 'Varnuc', 'Varcyt', ...
	'Anuc', 'Acyt', 'Anuc_ch', 'Acc','Bnuc', 'Bcyt', 'Bnuc_ch', 'Bcc',...
	'phinuc', 'phicyt', 'D2nuc','Q1','muQ1','mu1','mu2','mu3','sigma1','sigma2','sigma3',...
	'n2','n3','QbarNB','BbarNB','S_factor'};
nvars = length(varnames);

for j = 1:nvars
	if j <= 4 || isfield(Soln,varnames{j})
		eval([varnames{j},'_full = NaN(max(nt)+1,n_ncs,n_embryos);'])
	end
end

dD2nuc_full = NaN(max(nt)+1,n_ncs,n_embryos); % special case, 
% because fitting diffusion in the 2-cpt model is fraught, so we need to
% keep track of individual errorbars

for iii = 1:n_embryos

	t1 = Soln(iii).t;
	t14_0 = min(t1(:,end));

	%
	% Preallocate.
	%
	% The first four are special cases in which the variable name here does
	% not match the fieldname in Soln
	%
	Inuc1 = Soln(iii).nucsignal;
	Icyt1 = Soln(iii).cytsignal;
	Varnuc1 = Soln(iii).nucvar;
	Varcyt1 = Soln(iii).cytvar;
	for j = 5:nvars
		if exist([varnames{j},'_full'],'var')
			eval([varnames{j},'1 = Soln(iii).',varnames{j},';'])
		end
	end
	Anuc1(Anuc1 < 1e-5) = NaN;
	Acyt1(Acyt1 < 1e-5) = NaN;
	Anuc_ch1(Anuc_ch1 < 1e-5) = NaN;
	Acc1(Acc1 < 1e-5) = NaN;

	fitnuc = Soln(iii).fitnuc; % for the errorbars on D2

	%
	% Process the stages
	%
	for i = 1:n_ncs

		% -----------------------------------------------------------------
		% Normalizing the length of the nuclear cycles
		% -----------------------------------------------------------------
		
		%
		% Take the time points in the current NC and normalize them wrt
		% to the end of nc14 based on the average nc duration.
		%
		t_1 = t1(:,i); 
		if all(isnat(t_1)) % break out of NC's that have no data
			continue
		end
		tplot = minutes(t_1-t14_0);
		[tplot,iplot] = sort(tplot);
		vnan = isnan(tplot);
		tplot(vnan) = [];
		iplot(vnan) = [];

		if i < n_ncs
			m = (tendmean(i) - tstartmean(i))/(tend(iii,i) - tstart(iii,i));

		% special cases for nc 14, because not all gastrualtion captured:
		elseif ~isnan(tend(iii,4)) && ~isnan(tend(iii,3))
			% This case assumes there exist data in both nc 12 and nc
			% 13...this could be a problem because some embryos have only
			% nc 13,14, and others only nc14. In that case, we fall back on
			% alternatives.
			m = (tendmean(4) - tendmean(3))/(tend(iii,4) - tend(iii,3));

		elseif ~isnan(tend(iii,4))
			% The case for embryos that have nc13 but not 12:
			m = (tstartmean(5) - tstartmean(4))/(tstart(iii,5) - tstart(iii,4));

		else
			% final alternative: if there's only nc14...don't change
			% anything
			m = 1;			
		end
		t_hat = m*(tplot - tstart(iii,i)) + tstartmean(i);

		%
		% Time shift for nc 14
		%
		if i == n_ncs
			t00 = t_mesh(1,end);
		else
			t00 = 0;
		end
		if all(isnan(t_hat)) % break out of NC's that have no data
			continue
		end

		%
		% Calc the errorbars for D2 for this stage
		%
		nfiles = length(fitnuc{i});
		count_t = 1;
		dD2nuc1 = NaN(size(D2nuc1(:,i)));
		for j = 1:nfiles
			dD1 = fitnuc{i}{j}.errorbar68_2(3,:)';
			nt1 = length(dD1);
			dD2nuc1(count_t:nt1+(count_t-1)) = dD1;

			count_t = count_t + nt1 + 1; % The extra "+1" is to skip one
			% point to leave a NaN to separate time points, as we usually
			% do.
		end



		% -----------------------------------------------------------------
		% Do the interpolation. If only one data point, no interpolation
		% needed
		% -----------------------------------------------------------------
		if length(t_hat) > 1

			for j = 1:nvars
				if exist([varnames{j},'_full'],'var')
					eval([varnames{j},'_full(1:nt(i),i,iii) = ',...
						'interp1(t_hat+t00,',varnames{j},'1(iplot,i),t_mesh(1:nt(i),i));'])
				end
			end

			dD2nuc_full(1:nt(i),i,iii) = interp1(t_hat+t00,dD2nuc1(iplot),t_mesh(1:nt(i),i)); % dD2nuc1 only exists at this NC

		else % one data pt, no interpolation needed
			for j = 1:nvars
				if exist([varnames{j},'_full'],'var')
					eval([varnames{j},'_full(1:nt(i),i,iii) = ',...
						varnames{j},'1(iplot,i);'])
				end
			end
			dD2nuc_full(1:nt(i),i,iii) = dD2nuc1(iplot); % dD2nuc1 only exists at this NC

		end
	end
end


%
% Remove embryos that are fully NaN
%
v = all(squeeze(all(isnan(Anuc_full))))' | isnan(s);

for j = 1:nvars
	if exist([varnames{j},'_full'],'var')
		eval([varnames{j},'_full(:,:,v) = [];'])
	end
end
dD2nuc_full(:,:,v) = [];
s(v) = [];
n_embryos = size(Inuc_full,3);



%}





%% ========================================================================
% Average all variables and calc errorbars.
% =========================================================================
% {


if sem_flag
	D = sqrt(n_embryos);
else
	D = 1;
end

%
% Create average variables (but apparently only for Inuc?)
%
for j = 1:nvars
	if exist([varnames{j},'_full'],'var')
		if strcmp(varnames{j},'Inuc')
			eval([varnames{j},' = mean(',varnames{j},'_full,3,''omitmissing'');'])
			eval([varnames{j},'_e = std(',varnames{j},'_full,[],3,''omitmissing'')/D;'])
		end
	end
end

% special case because really need weighted mean for the diffusivity
w = 1./dD2nuc_full./sum(1./dD2nuc_full,3,'omitmissing');
D2nuc = sum(w.*D2nuc_full,3,'omitmissing');

std_weighted = w.*(D2nuc_full - D2nuc).^2;
D2nuc_e = sqrt(sum(std_weighted,3,'omitmissing')./sum(w,3,'omitmissing'))/D;



%}



%% ========================================================================
% Create derived quantities
% =========================================================================
% {



%
% Major derived quantities
%
Cnuc_full = 1./(V_eff*Anuc_full)/ 0.602;
Ccyt_full = 1./(V_eff*Acyt_full)/ 0.602;
Cnuc_ch_full = 1./(Vnuc_ch*Anuc_ch_full)/ 0.602;

phicc_full = Vcc./Vnuc_ch.*Acc_full./Anuc_ch_full;
phiuncorrel_full = phinuc_full-phicc_full;

Cfree_full = Cnuc_full.*(1-phinuc_full);
CDNAbound_full = Cnuc_full.*phicc_full;
Cuncorrel_full = Cnuc_full.*(phinuc_full - phicc_full);

%}



%% ========================================================================
% Bin quantities into AP locations, averaged in mid nc14
% =========================================================================
% {
%
%
%

%
% Create bin size if not specified, then create mesh
%
s1 = min(s);
ds = 0.06;
s_mesh = (s1:ds:max(s))';
sb1 = max(0,s_mesh - ds/2);
sb2 = s_mesh + ds/2; sb2(end) = max(s)*(1+1e-4);
ns = length(sb1);

%
% Pre-allocate
%
varnames2 = {'Cnuc','Ccyt','Cnuc_ch','phicc','phiuncorrel','Cfree','CDNAbound','Cuncorrel',...
	'NCRI'};
Varnames = [varnames varnames2]; % combine variable names
for j = 1:length(Varnames)
	if exist([Varnames{j},'_full'],'var')
		eval([Varnames{j},' = zeros(ns,1);'])
		eval([Varnames{j},'_e = zeros(ns,1);'])
	end
end


%
% Do averaging in time (nc14 from 30 min to 10 min bf gast) and binning in
% space (along s_mesh)
%
vt = t(:,NC_list == 14)+t0 > -30 & t(:,NC_list == 14)+t0 < -10;
n = zeros(ns,1);
for i = 1:ns
	v = s >= sb1(i) & s < sb2(i);
	if sum(v) <= 1
		v = false(size(v));
	end
	n(i) = sum(v);

	if sem_flag
		D = sqrt(n(i));
	else
		D = 1;
	end
	% D = sqrt(n(i)); % assume we're doing SEM


	for j = 1:length(Varnames)
		if exist([Varnames{j},'_full'],'var')
			% extract the relevant time points of embryos that are of the
			% right orientation:
			x = squeeze(eval([Varnames{j},'_full(vt,NC_list == 14,v)'])); 

			% correction if there are no valid points during mid nc14
			vnan = all(isnan(x));
			for k = 1:length(vnan)
				if vnan(k)
					vnan2 = find(v);
					x1 = squeeze(eval([Varnames{j},'_full(:,NC_list == 14,vnan2(k))']));
					n1 = min(length(x1),size(x,1));
					x(1:n1,k) = x1(1:n1);
				end
			end
			
			eval([Varnames{j},'(i) = mean(x(:),''omitmissing'');'])
			eval([Varnames{j},'_e(i) = std(x(:),''omitmissing'')/D;'])
		end
	end
end



%}



%% ========================================================================
% Save output variables
% =========================================================================
% {

%
% Header variables
%
data.t = t;
data.s = s_mesh;
data.s_list = s;
data.n = n;
data.t0 = t0;
data.idx_de = idx_de;
data.V_eff = V_eff;
data.Vcc = Vcc;
data.Vnuc_ch = Vnuc_ch;

%
% Averaged variables
%
varnames2 = {'Cnuc','Ccyt','Cnuc_ch','phicc','phiuncorrel','Cfree','CDNAbound','Cuncorrel',...
	'NCRI'};
Varnames = [varnames varnames2]; % combine variable names
for j = 1:length(Varnames)
	if exist([Varnames{j},'_full'],'var')
		data.(Varnames{j}) = eval(Varnames{j});
	end
end


%
% errorbars
%
for j = 1:length(Varnames)
	if exist([Varnames{j},'_full'],'var')
		data.([Varnames{j},'_e']) = eval([Varnames{j},'_e']);
	end
end


%
% Full variables (all embryos), just for the "basic" variables + Cnuc
%
for j = 1:nvars
	if exist([varnames{j},'_full'],'var')
		data.allembryos.(varnames{j}) = eval([varnames{j},'_full']);
	end
end
data.allembryos.Cnuc = Cnuc_full;



%}





%% ========================================================================
% Plot orientation, if asked for
% =========================================================================
% {


if ~exist('yesplot','var') || isempty(yesplot) || ~yesplot || yesplot == 0
	return
end

%
% Do the plotting
%
figure('pos',[9         102        1899         592])
TL = tiledlayout(2,5,'TileSpacing','Compact','Padding','Compact');
C = colormap('lines');



%
% Make cell variables for each of the plots (X,Y,S)
%
Y = {Cnuc, Cfree, CDNAbound, Cuncorrel, CDNAbound
	Cnuc_ch, 1-phinuc, phicc, phinuc-phicc, Cuncorrel};

Sx = {0, 0, 0, 0, Cfree_e
	0, 0, 0, 0, Cfree_e};
Sy = {Cnuc_e, Cfree_e, CDNAbound_e, Cuncorrel_e, CDNAbound_e
	Cnuc_ch_e, phinuc_e, phicc_e, sqrt(phinuc_e.^2+phicc_e.^2), Cuncorrel_e};


%
% Make graph-labeling variables
%
Xlabel = [repmat({'AP coord'},1,4),{'Free nuc conc'}
	repmat({'AP coord'},1,4),{'Free nuc conc'}];
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
		errorbar(s_mesh,Y{i,j},Sy{i,j},'o-','linewidth',2)
		ylim([0 YLIM2(i,j)])
	else

		errorbar(Cfree,Y{i,j},Sy{i,j},Sy{i,j},Sx{i,j},Sx{i,j},'o-','linewidth',2)
		
	end
	xlabel(Xlabel{i,j})
	ylabel(Ylabel{i,j})
	set(gca,'fontsize',12)

	end
end





%}

