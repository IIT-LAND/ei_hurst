function B_invivo_1_DREADDpfc_excitation

% add nonfractal toolbox to MATLAB path
addpath /Users/mlombardo/Dropbox/matlab/nonfractal-master/m;

%% load in pheno data
rootpath = '/Users/mlombardo/Dropbox/Manuscripts/AIMS_Hurst_Sex/reproAnalysis/DREADDexcitation';
datapath = '/Users/mlombardo/Dropbox/data/ag_EI_DREADD_data/DREADDexcitiation/raw';
pheno_path = fullfile(rootpath,'pheno');
pdata = readtable(fullfile(pheno_path,'tidy_pheno_data.csv'),'delimiter',',');


%% Load in data
sublist = pdata.filename;
ntimepoints = 4040;
tidy_data = zeros(length(sublist),ntimepoints);
for i = 1:length(sublist)
    fname = fullfile(datapath,sublist{i});
    tmp_data = readtable(fname,'ReadVariableNames',false);
    tidy_data(i,:) = tmp_data.Var1';
end

for i = 1:ntimepoints
    volnames{i} = sprintf('vol%03d',i);
end
tidydata2write = cell2table(num2cell(tidy_data),'VariableNames',volnames);
writetable(tidydata2write,fullfile(rootpath,'tidy','tidy_data.csv'),'FileType','text','delimiter',',');


%% compute H on tidy_data_long
Hfilter = 'haar';
waveletType = 'dwt';
lb_param = [-0.5,0];
ub_param = [1.5,10];

window_length = 512;
start_idx = 1:((ntimepoints)-(window_length-1));
end_idx = 512:(ntimepoints);

% compute H
Hwin = zeros(size(tidy_data,1),length(start_idx));

parfor iwindow = 1:length(start_idx)
    disp(iwindow);
    idx2use = start_idx(iwindow):end_idx(iwindow);
    data2use = tidy_data(:,idx2use);
    Hwin(:,iwindow) = bfn_mfin_ml(data2use', ...
        'filter',Hfilter, ...
        'lb',lb_param,'ub',ub_param, ...
        'wavelet',waveletType, ...
        'verbose',false);
end % for iwindow

%% compute fALFF on tidy_data_long

% compute fALFF
fALFFwin = zeros(size(tidy_data,1),length(start_idx));

for isub = 1:length(sublist)
    for iwindow = 1:length(start_idx)
        idx2use = start_idx(iwindow):end_idx(iwindow);
        data2use = tidy_data(isub,idx2use);
        fALFFwin(isub,iwindow) = falff(data2use);
    end % for iwindow
end % for isub

%%
for i = 1:size(Hwin,2)
    volnames3{i} = sprintf('window_%04d',i);
end

Hwin_tab = cell2table(num2cell(Hwin),'VariableNames',volnames3);
pdata_final = [pdata Hwin_tab];

fname2write = fullfile(rootpath,'pheno','pheno_data+Hwin_dreaddexcitation.csv');
writetable(pdata_final,fname2write,'FileType','text','delimiter',',');


fALFFwin_tab = cell2table(num2cell(fALFFwin),'VariableNames',volnames3);
pdata_final = [pdata fALFFwin_tab];

fname2write = fullfile(rootpath,'pheno','pheno_data+fALFFwin_dreaddexcitation.csv');
writetable(pdata_final,fname2write,'FileType','text','delimiter',',');
end % function B_invivo_1_DREADDpfc_excitation


%%
function result = falff(data)
%
%   data = 1:ntimepoints vector
%

%%
TR = 1;
Fs = 1/TR; % sampling rate in Hz
bp = [0.01,0.03]; % infra-slow band
fp = [0.01, 0.1]; % full band

Nobs = size(data,2); % number of timepoints

% mean center
mean_data = nanmean(data);
mc_data = data - mean_data;

% run fft
xdft = fft(mc_data);
xdft = xdft(1:Nobs/2+1);
psdd1(1,:) = (1/(Fs*Nobs)) * abs(xdft).^2;
psdd1(1,2:end-1) = 2*psdd1(1,2:end-1);
freqs1 = 0:Fs/Nobs:Fs/2;

pband1 = bandpower(psdd1,freqs1,bp,'psd');
ptot1 = bandpower(psdd1,freqs1,fp,'psd');

% falff
result = pband1./ptot1;
end % function falff
