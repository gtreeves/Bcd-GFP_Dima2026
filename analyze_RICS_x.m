function Gs = analyze_RICS_x(filename,data_ch,genotype,w0,wz,varargin)
%Calculates the RICS experimental ACF for total, nuc, and cyt.
%
%function Gs = analyze_RICS_x(filename,data_ch,genotype,w0,wz,varargin)
%
% This function is "x" because we are dividing the image into multiple
% spatial domains. We assume the image is rectangular with an aspect ratio
% of roughly 4, so we are dividing it into 4 near-squares.
%
% This function first reads in an lsm-file of a time series acquisition
% suitable for RICS analysis.  Next, after image ROI is cropped and the
% image is background subtracted, the ACF is calculated from the data.  (A
% mask may be calculated first, if asked for).
%
% Inputs:
% "filename": the name of the file, including the path to it
% "data_ch": the channel of your image stack where the data (on which RICS
%	should be performed) reside. Alt: the wavelength of the laser that
%	excited the data channel.
% "genotype": a descriptive genotype of the biological sample
% "w0": the e^(-2) xy waist of the Gaussian PSF
% "wz": the axial waist 
%
% Optional argument varargin can consist of these things, in this order:
%	* "sws": input that tells the function whether to use a sliding window to
%		subtract off an average frame, and/or how many frames go into the
%		sliding window. This input can be logical or numeric. If this input
%		is a logical, it will dictate whether to subtract off the average
%		frame.  By "average frame" we mean the pixel-by-pixel average of
%		all frames in a given window around frame "i". You average the
%		frames in this window in the t-direction and you get an H-by-W
%		average image.  This is subtracted off of frame "i" IM. Then, the
%		average intensity of that average frame is added back so that the
%		value of your function will not be very close to zero. This value
%		of this procedure is it will get rid of the effect of large-scale
%		structures in your image, such as the presence of a nucleus that
%		has a different brightness than the cytoplasm.  This is according
%		to Digman et al., 2005, Biophys J, 88: L33-L36. If this is given as
%		logical true, then by default the function will subtract off a
%		window of 5 frames on either side of frame "i". If given as
%		numeric, it will be the number of frames on either side of frame
%		"i" that you'd like to include in a moving average. If you specify
%		this number to be negative, then the average of the entire time
%		series stack is subtracted off. If a number greater than zero, then
%		an average frame in a sliding window centered at point "i" will be
%		subtracted from frame "i".  The total number of frames in the
%		sliding window average is 2*(sliding_window_size) + 1. Default, 5.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "mask_ch": numeric input that tells the function which channel has
%		the image data from which to derive the mask.  This can be used to
%		calculate the autocorrelation function for any arbitrarily-shaped
%		region, and perhaps even disjoint regions.  If passed as a scalar
%		numerical zero, or logical false, then the mask will not be used
%		and the autocorrelation function for the entire frame will be
%		calculated. If passed as a low-value integer (like 1,2,3,4, etc),
%		then that will be the numerical index of the mask channel. If
%		passed as a larger number, it will specify the wavelength in nm for
%		the laser that was used to image the mask channel. Default,
%		whatever channel corresponds to a laser wavelength in the 500s, if
%		such a channel exists. If not, then this argument defaults to
%		false.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "usemask": input that tells the function whether to use a mask as an
%		ROI to calculate the data-ACF.  This input can be logical or
%		numeric. If this input is a logical, it will dictate whether a mask
%		is used at all. If so, then by default 10 frames are averaged
%		together to get a more robust mask. If it is numeric, it dictates
%		the number of frames you'd like to average together to better
%		calculate the mask. If a mask is calculated then both a cytoplasmic
%		and nuclear mask will be calculated and used to get the nuclear and
%		cytoplasmic intensities. Note that this is contingent on there
%		being a "mask_ch". Default, 10.
%
%		NOTE: This argument can also be a group index vector that says we want
%		the windowed average frame subtraction to be the exact same frames as
%		the ones grouped together for usemask.
%
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "zerospad": Whether or not you'd like to choose to pad your frame
%		with zeros so that the correlation does not use periodic
%		wrap-around (thereby resulting in non-physical terms in your
%		correlation function).  While having non-physical terms sounds
%		scary, it actually does not seem to affect much.  Default, false.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "hasbackground": Logical variable that tells whether the edge of the
%		embryo is present in the image, and if so, that we need to
%		calculate the border between embryo and background (outside the
%		embryo). Default, false.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "xilimit": Integer that describes how much of the ACFs and CCF we
%		want to save, since most of it is just background levels (near zero).
%		Default, 51.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "yesplot": whether you want to plot the outcome.  Default, "false".  
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%
%
% Outputs:
% The function returns the structures 'Gs' and 'Gs_t with fields specified
% in the "Middle remarks" and at the end of the file. In, "Gs" is the
% standard structure file containing the ACFs/CCF of interest (all time
% stamps/groups averaged together), while "Gs_t" keeps the groups separated
% so that a time course can be constructed.



%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	sws = varargin{iArg}; else
	sws = 5;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	mask_ch = varargin{iArg}; else 
	mask_ch = [];
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	usemask = varargin{iArg}; else 
	usemask = 10;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	zerospad = varargin{iArg}; else
	zerospad = false;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	hasbackground = varargin{iArg}; else 
	hasbackground = false;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
 	xilimit = varargin{iArg}; else
 	xilimit = 51;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
 	yesplot = varargin{iArg}; else
 	yesplot = false;
end%, iArg = iArg + 1;


%
% filename and load in image and metadata
%
if ~ischar(filename)
	error('wrong argument type')
end
v = find(filename == '\' | filename == '/');
v = v(end);
filenameshort = filename(v+1:end);
pth = filename(1:v);

stabilize = 150; % minimum frame size allowed
bg1 = true; % subtract background within ftn_imread (calls ftn_calcbg)
[IM,lsminf1,lsminf2,bg_lvl,sig0,istart,iend,Grpidx,tstart] = ftn_imread(filename,bg1,stabilize);
if Grpidx
	usemask = Grpidx;
end

[H,W,num_ch,D,num_frames] = size(IM);
% We are currently assuming D = 1, so we need to squeeze IM.
IM = squeeze(IM);

%
% Mask ch and data ch override. This works only for the czi files, because
% I don't know if I can trust the metadata of the lsm files. Luckily, all
% of our lsm files have data_ch = 1.
%

% define lambda
if isfield(lsminf2,'cziinfo')
	cziinf = lsminf2.cziinfo.cziinf.Channel;
	fnames = fieldnames(cziinf);
	v1 = find(contains(fnames,'ExcitationWavelength'));

	if ~isempty(v1)
		lambda = zeros(length(v1),1);
		for i = 1:length(v1)
			lambda(i) = cziinf.(fnames{v1(i)});

		end
	end
else
	lambda = [488 590];
end

if isempty(data_ch)
	data_ch = 488;
end
if isempty(mask_ch)
	mask_ch = 590;
end

% define data channel
if isnumeric(data_ch) && data_ch > 10
	data_ch = find(lambda > data_ch-35 & lambda < data_ch+35);

elseif ~isnumeric(data_ch)
	data_ch = 1; % default
end

% define mask channel
if isnumeric(mask_ch) && mask_ch > 10
	mask_ch = find(lambda > mask_ch-35 & lambda < mask_ch+35);

elseif ~isnumeric(mask_ch)
	mask_ch = 2; % default
end


%
% Check for the "channels" input
%
if isempty(data_ch) || data_ch > num_ch || data_ch < 1
	fprintf('Number of channels in the z-stack: %d\n',num_ch)
	error(['Make sure your "data_ch" input is included between 1 and ',num2str(num_ch)])
end
if (length(usemask) > 1 || usemask) && (isempty(mask_ch) || mask_ch > num_ch || mask_ch < 1)
	fprintf('Number of channels in the z-stack: %d\n',num_ch)
	error(['Make sure your "mask_ch" input is included between 1 and ',num2str(num_ch)])
end

%
% Size and scalings
%
scalings = lsminf2.VOXELSIZES*1e6; % microns/pixel
dr = scalings(1);

%
% Determine the orientation of the rectangle and the number of groupings in x
%
if H < W
	iswide = true;
else
	iswide = false;
end
nx = max(round(W/H),round(H/W)); % number of groupings in x
if nx == 1
	error('Image aspect ratio close to one. Exiting program')
end

%
% Pixel, line, and frame time
%
taup = lsminf2.ScanInfo.PIXEL_TIME{data_ch}*1e-6; % pixel dwell time in s
H_orig = lsminf2.DimensionY;
tauL = lsminf2.TimeInterval/H_orig; % line time in s
tauf = lsminf2.TimeInterval; % total frame time in s


%
% Middle remarks:
%
Gs.pth = pth; % path to image file
Gs.filenameshort = filenameshort; % just the file name
Gs.filename = filename; % path + filename
Gs.genotype = genotype; 
Gs.side = ''; % placeholder
Gs.metadata.data_ch = data_ch;
Gs.metadata.mask_ch = mask_ch;
Gs.metadata.istart = istart; % start index (or [scene,index])
Gs.metadata.iend = iend; % end index (or [scene,index])
Gs.metadata.lsminf1 = lsminf1; % empty; deprecated
Gs.metadata.lsminf2 = lsminf2; % structure with more detailed metadata
Gs.metadata.scalings = scalings; % pixels per micron
Gs.metadata.bg = bg_lvl; % the background intensity in absence of signal
Gs.metadata.std_bg = sig0; % std in intensity due to dark current
Gs.metadata.stabilize = stabilize; % flag indicating whether image was stabilized upon read
Gs.metadata.S_factor = []; % relationship between digital level and photons
Gs.metadata.subtravg = sws > 0 || sws == -1; % flag stating if sliding window subtraction was used
Gs.metadata.sliding_window_size = sws; % side of sliding window used for subtraction
Gs.metadata.zerospad = zerospad; % flag stating whether zeros were used to pad the image before FFT, nearly always false
Gs.metadata.usemask = usemask; % flag stating whether nuclear mask was used
Gs.metadata.hasbackground = hasbackground; % flag stating whether there is any non-embryo part of the image, nearly always false
Gs.metadata.xbg = NaN; % placeholder
Gs.metadata.ybg = NaN; % placeholder
Gs.metadata.H = H; % number of rows in image (after stabilization, if used)
Gs.metadata.W = W;  % number of cols in image (after stabilization, if used)
Gs.metadata.z_depth = D; % number of z-slices in image (usually 1 for RICS)
Gs.metadata.num_frames = num_frames; % number of frames in time course
Gs.metadata.num_ch = num_ch; % number of color channels
Gs.metadata.dr = dr; % pixel size (microns per pixel)
Gs.metadata.w0 = w0; % PSF waist size in xy plane
Gs.metadata.wz = wz; % PSF waist size along axial (z) direction
Gs.metadata.taup = taup*1e6; % pixel time in microseconds
Gs.metadata.tauL = tauL*1000; % line time in ms
Gs.metadata.tauf = tauf; % frame time in s




% =========================================================================
% Now that the pre-processing is done, it is time to do the analysis here
% in the second half of this function.
% =========================================================================


% 
% Find the mask, if asked for.  Note that this part is not well optimized.
% 
if (length(usemask) > 1 || usemask > 0) && isnumeric(mask_ch) && mask_ch > 0 && mask_ch <= num_ch
	if hasbackground
		[xbg,ybg,Vbg] = find_edge(IM,Gs,0);
		Gs.metadata.xbg = xbg;
		Gs.metadata.ybg = ybg;
	else		
		Vbg = [];
	end
	if isnumeric(usemask)
		% then usemask will be either the number of slices we want to group
		% together (nslices), or it could be a vector of group indices.
		[nucmask,cytmask,~,~,Grpidx] = find_mask(IM,Gs,Vbg,0,usemask,true);
		usemask = true;
	else
		[nucmask,cytmask,~,~,Grpidx] = find_mask(IM,Gs,Vbg,0,[],true);
	end
else
	nucmask = NaN;
	cytmask = NaN;
end



%
% Preallocate in advance of the coming for loop
%
ACFnames = {'Gstot','Gsnuc','Gscyt','Gscc','Gsnuc_ch'}; % xilimit-by-xilimit-by-nx

scalarnames = {'totsignal','totvar','nucsignal','nucvar',...
	'cytsignal','cytvar','nuc_chsignal','nuc_chvar'}; % nx-by-1

for k = 1:length(ACFnames)
	eval([ACFnames{k},' = NaN(xilimit,xilimit,nx);'])
end
for k = 1:length(scalarnames)
	eval([scalarnames{k},' = NaN(nx,1);'])
end

%
% Establish a for loop that will separately process different spatial
% grides within the image. The image should be a rectangle with the long
% dimension close to an integer multiple of the short direction, so we will
% process the image into a series of near-squares.
%
for j = 1:nx
	if iswide
		x_idx = round((j-1)*W/nx+1):round(j*W/nx);
		I = double(squeeze(IM(:,x_idx,data_ch,:)));
		I2 = double(permute(IM(:,x_idx,:,:),[1 2 4 3]));
		I_nuc = double(squeeze(IM(:,x_idx,mask_ch,:)));
		nucmask1 = nucmask(:,x_idx,:);
		cytmask1 = cytmask(:,x_idx,:);
	else
		x_idx = round((j-1)*H/nx+1):round(j*H/nx);
		I = double(squeeze(IM(x_idx,:,data_ch,:)));
		I2 = double(permute(IM(x_idx,:,:,:),[1 2 4 3]));
		I_nuc = double(squeeze(IM(x_idx,:,mask_ch,:)));
		nucmask1 = nucmask(x_idx,:,:);
		cytmask1 = cytmask(x_idx,:,:);
	end

	%
	% Calculate the autocorrelation function from our data
	%
	if sws || sws > 0
		subtravg = true;
	end
	[G,Imean_tot,Var_tot] = make_ACF_data(I,false,subtravg,sws,Grpidx);
	Gstot(:,:,j) = G(1:xilimit,1:xilimit);
	totsignal(j) = mean(Imean_tot); totvar(j) = mean(Var_tot);
	if usemask
		[G,Imean_nuc,Var_nuc] = make_ACF_data(I,nucmask1,subtravg,sws,Grpidx);
		Gsnuc(:,:,j) = G(1:xilimit,1:xilimit);
		nucsignal(j) = mean(Imean_nuc); nucvar(j) = mean(Var_nuc);

		[G,Imean_cyt,Var_cyt] = make_ACF_data(I,cytmask1,subtravg,sws,Grpidx);
		Gscyt(:,:,j) = G(1:xilimit,1:xilimit);
		cytsignal(j) = mean(Imean_cyt); cytvar(j) = mean(Var_cyt);
		if num_ch > 1 && mask_ch ~= data_ch

			%
			% Cross correlation
			%
			G = make_ACF_data(I2,nucmask1,subtravg,sws,Grpidx);
			Gscc(:,:,j) = G(1:xilimit,1:xilimit);

			%
			% Nuclear channel only
			%
			[G,Imean_nuc_ch,Var_nuc_ch] = make_ACF_data(I_nuc,nucmask1,subtravg,sws,Grpidx);
			Gsnuc_ch(:,:,j) = G(1:xilimit,1:xilimit);
			nuc_chsignal(j) = mean(Imean_nuc_ch); nuc_chvar(j) = mean(Var_nuc_ch);
		end
	end
end

%
% Average fluorescence intensity and variance of each pixel in the time
% course (for N+B)
%
I0 = double(squeeze(IM(:,:,data_ch,:)));
I = mean(I0,3);
V = var(I0,[],3);


%
% Closing remarks
%
Gs.Gstot = Gstot;
Gs.Gsnuc = Gsnuc;
Gs.Gscyt = Gscyt;
Gs.Gscc = Gscc;
Gs.Gsnuc_ch = Gsnuc_ch;
Gs.nucmask = nucmask;
Gs.cytmask = cytmask;
Gs.totsignal = totsignal;
Gs.nucsignal = nucsignal;
Gs.cytsignal = cytsignal;
Gs.nuc_chsignal = nuc_chsignal;
Gs.totvar = totvar;
Gs.nucvar = nucvar;
Gs.cytvar = cytvar;
Gs.nuc_chvar = nuc_chvar;
Gs.Intensity = I;
Gs.Variance = V;

% Optional to auto-save to mat file. Could take up a lot of space.
% save([Gs.filename(1:end-4),'_Gs.mat'],'Gs');
% save([Gs.filename(1:end-4),'_Gs_t.mat'],'Gs_t');











