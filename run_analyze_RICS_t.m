function [Soln,Gst_out] = run_analyze_RICS_t(filenames,stages,istarts,iends,...
	embryo_indices,basenames,data_ch,genotype,w0,wz,varargin)
%Runs "analyze_RICS" on each czi file as specified in "filenames"
%
%function Soln = run_analyze_RICS_t(filenames,stages,istarts,iends,...
%	 embryo_indices,basenames,data_ch,genotype,w0,wz,varargin)
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
%		calculate the mask. If a mask is calculated, then both a cytoplasmic 
%		and nuclear mask will be calculated and used to get the nuclear and
%		cytoplasmic intensities. Note that this is contingent on there
%		being a "mask_ch". Default, 10.
%
%		NOTE: in the future, I want to have an argument that says we want
%		the sws to be the exact same frames as the ones grouped together
%		for usemask. That is tricky because the sws is a sliding window,
%		while the frames grouped for usemask are static (not sliding).
%
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "outfilename": If you want the output save files to have a name other
%		than just "run_analyze_xsCurrent.mat", it appends this
%		character variable onto the end of the name.  Must be a character
%		variable, otherwise defaults.  Default, empty string, ''.
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
%	* "yesplot": whether you want to plot the outcome.  Default, "false".
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.


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
	outfilename = varargin{iArg}; else
	outfilename = '';
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
	yesplot = varargin{iArg}; else
	yesplot = false;
end%, iArg = iArg + 1;


if ~ischar(outfilename)
	outfilename = '';
end
if ~isempty(outfilename) && ~strcmp(outfilename(end),'_')
	outfilename(end+1) = '_';
end

NC_list = 10:14;
clk1 = clock; %#ok<CLOCK>
clk = [num2str(clk1(1)),'-',num2strDU(clk1(2),2),'-',num2strDU(clk1(3),2),'_',...
	num2strDU(clk1(4),2),'-',num2strDU(clk1(5),2),'-',num2strDU(round(clk1(6)),2)];

%
% Fitting params
%
B0 = [0 -1e-3 1e-3]; % set B0 to zero so we don't overfit if it drops below zero
D20 = 3; % just for the 2-cpt model for now
W0 = []; % allow w0 to vary slightly for robustness of fit


%
% Loop over each embryo index
%
Soln = [];
Gst_out = [];
n_embryos = max(embryo_indices);

for iii = 1:n_embryos % but doing multiple embryos is not recommended in case of crash and memory issues
	k_embryo = find(embryo_indices == iii);
	n_files = length(k_embryo);

	%
	% Loop over every filename/stage for this embryo
	%
	GS_t = cell(1,length(NC_list));
	Data = cell(1,length(NC_list));
	tstamp = cell(1,length(NC_list));
	for ii = 1:n_files

		%
		% Filename, stage, istart, iend for this filename/stage
		%
		filename1 = filenames{k_embryo(ii)};
		stage1 = stages(k_embryo(ii));
		istart1 = istarts(k_embryo(ii),:);
		iend1 = iends(k_embryo(ii),:);
		basename1 = basenames{k_embryo(ii)};

		%
		% Run the image analysis
		%
		[~,Gs_t] = analyze_RICS_t(filename1,data_ch,genotype,w0,wz,...
			sws,mask_ch,usemask,istart1,iend1,basename1,...
			zerospad,hasbackground,yesplot);
		data1 = fit_both_models_2step(Gs_t,B0,D20,W0);
		tstamp1 = data1.t(1);

		%
		% Place our analysis into the appropriate stage bin
		%
		vstage = find(stage1 == NC_list);
		if isempty(Data{vstage})
			GS_t{vstage} = Gs_t;
			Data{vstage} = data1;
			tstamp{vstage} = tstamp1;
		else
			Gs_t2 = [GS_t{vstage} Gs_t];
			data2 = [Data{vstage} data1];
			tstamp2 = [tstamp{vstage} tstamp1];
			[tstamp2,isort] = sort(tstamp2);
			tstamp{vstage} = tstamp2;
			Data{vstage} = data2(isort);
			GS_t{vstage} = Gs_t2(isort);
		end


		disp(['ii = ',num2str(ii),' out of ',num2str(n_files),' files/stages'])

	end

	%
	% For this embryo, combine all structures from "Data" into one
	% megastructure!
	%
	data = timecourse_structurecombine(Data,NC_list);
	Soln = [Soln;data]; %#ok<AGROW>
	Gs_t = timecourse_structurecombine(GS_t,NC_list);
	Gst_out = [Gst_out;Gs_t]; %#ok<AGROW>

	if ~exist('Mat','dir')
		mkdir('Mat')
	end

	if iii == n_embryos
		matfilename1 = ['Mat/',clk,'_',outfilename,'timecourse_Soln'];
		matfilename2 = ['Mat/',clk,'_',outfilename,'timecourse_Gst'];
	else
		matfilename1 = ['Mat/',clk,'_',outfilename,'timecourse_Soln_temp'];
		matfilename2 = ['Mat/',clk,'_',outfilename,'timecourse_Gst_temp'];

	end
	save(matfilename1,'Soln')
	save(matfilename2,'Soln','Gst_out','-v7.3')
	disp(['iii = ',num2str(iii),' out of ',num2str(n_embryos),' embryos'])

	%
	% Optional: save time course file to original path. Takes up a lot of
	% space
	%
	% kslash = find(filename1 == '/' | filename1 == '\');
	% pth1 = filename1(1:kslash(end));
	% if isempty(basename1), u = ''; else, u = '_'; end
	% save([pth1,basename1,u,'fulltimecourse'],'data','Gs_t','-v7.3');
end











