function [filenames,stages,istarts,iends,embryo_indices,basenames,pths_u,pths,filenamesshort] = ...
	timecourse_parser(filename,sheetname)
%Examines the excel or csv file in "filename" and parses output
%
%function [filenames,stages,istarts,iends,embryo_indices,pths_u,pths,filenamesshort] = ...
%	timecourse_parser(filename,sheetname)
%
% The purpose of this function is to extract information from the excel
% file that describes the folders and filenames of our czi files, as well
% as the nuclear cycle and start and end time points. The output of this
% file is a list of different files to run through analyze_RICS_t, as well as
% the start and end frame for each nuclear cycle in each file. This
% function can also read-in a csv file, which could be stored in a given
% folder. So an alternative input, besides an excel filename, is either a
% csv file, or even just a folder. If it's just a folder, this file will
% scan the folder for an appropriate excel, text, or csv file to parse.
%
% Optional argument "sheetname": Character variable of the sheetname that
% you're looking into in the given excel file. If passed as empty, by
% default, the first sheet in the excel file will be used.


%
% Read data
%
opts = detectImportOptions(filename);
opts.DataRange = 'A2';
opts.VariableNamesRange = 'A1';
opts.Sheet = sheetname;
opts.VariableNamingRule = 'preserve';

try
	T = readtable(filename,opts);
	% T = readtable(filename,'Sheet',sheetname,'VariableNamingRule','preserve');
catch
	T = readtable(filename,opts,'VariableNamingRule','preserve');
end
fnames = fieldnames(T);

%
% Extract stages, start and end indices
%
istages = find(contains(fnames,'Stage') | strcmpi(fnames,'nc') | strcmpi(fnames,'mitosis'));
stages = T.(fnames{istages});

idx = istages+1:length(fnames);
V = NaN(length(stages),length(idx));
for i = 1:length(idx)
	a = T.(fnames{idx(i)});
	if isnumeric(a)
		V(:,i) = a;
	else
		break
	end
end
V(:,all(isnan(V))) = [];
n = size(V,2); % standard: n = 2, which is [istart,iend]. However, if n = 4,
% then maybe we have [start scene, start frame, end scene, end frame],
% which is what we'd get with a "RICS_tiles" timecourse. 
if n >= 4
	istarts = V(:,1:2);
	iends = V(:,3:4);
else % treat as if n == 2
	istarts = V(:,1);
	iends = V(:,2);
end

%
% Parse path and filename variables. As you go down each array, we will
% fill in any missing values with the value above it.
% 
P = table2array(T(:,1:istages-1));
for j = 1:size(P,2)
	for i = 2:size(P,1)
		if isempty(P{i,j})
			P{i,j} = P{i-1,j};
		end
	end
end

% double-checking for valid ending to each filename
for i = 1:size(P,1)
	if length(P{i,end}) <= 4 || (~strcmp(P{i,end}(end-4:end),'.czi') && ~strcmp(P{i,end}(end-4:end),'.lsm') && ...
			~strcmp(P{i,end}(end-4:end),'.lif'))
			P{i,end}(end+1:end+4) = '.czi'; % assume czi file
	end
end
filenamesshort = P(:,end);

% Concatenate paths, double-checking for a valid filesep
pths = cell(size(P,1),1);
filenames = cell(size(P,1),1);
for i = 1:size(P,1)
	for j = 1:size(P,2)-1
		if ~isempty(P{i,j}) && ~strcmp(P{i,j}(end),'/') && ~strcmp(P{i,j}(end),'\')
			P{i,j}(end+1) = filesep;
		end
	end
	pths{i} = horzcat(P{i,1:end-1});
	filenames{i} = horzcat(P{i,1:end});
end

[pths_u,~,ic] = unique(pths);


%
% Check to see if there is more than one embryo in each path. These are
% designated by "embryoX" or "embryoXX" or "embryo X" or "embryo XX"
%
embryo_indices = zeros(length(pths),1);
current_index = 1;

basenames = cell(length(pths),1);
for i = 1:length(pths_u)
	k = find(ic == i);
	
	basename = cell(length(k),1);
	for j = 1:length(k)
		filenamesshort1 = filenamesshort{k(j)};

		if contains(filenamesshort1,'embryo')
			if isspace(filenames{k(j)}(7)) % the first six characters are "embryo"
				if ~isnan(str2double(filenamesshort{k(j)}(8:9))) % then "embryo XX"
					k_end = 9;
				else % then "embryo X"
					k_end = 8;
				end
			else
				if ~isnan(str2double(filenamesshort{k(j)}(7:8))) % then "embryoXX"
					k_end = 8;
				else % then "embryoX"
					k_end = 7;
				end
			end

			basename{j} = filenamesshort{k(j)}(1:k_end); % something like "embryoX"
		else
			basename{j} = '';
		end
	end

	[b,~,ic2] = unique(basename);
	for j = 1:length(b) % looping over every base name in pths_u{i}
		embryo_indices(k(ic2 == j)) = current_index;
		basenames(k(ic2 == j)) = b(j);
		current_index = current_index + 1;
	end



end











