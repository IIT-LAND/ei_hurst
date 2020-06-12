function result = A_insilico_3_boldsim_neuralts_oof(MAKE_PLOT, fstem)
%
% INPUT
%   MAKE_PLOT = set to 1 to make a plot, otherwise set to 0.
%   fstem = 'oof';
%

addpath /Users/mlombardo/Dropbox/matlab/spm12;

%%
% BOLD hrf parameters
sampling_rate = 500;
dt = 1/sampling_rate;
p = [6,16,1,1,6,0,32];
[hrf, p] = spm_hrf(dt,p);

fmri_tp = 1800;

% parameters for H computation
lb_param = [-0.5,0];
ub_param = [1.5,10];
Hfilter_type = 'haar';

%%
datapath = '/Users/mvlombardo/projects/EI_hurst/reproAnalysis/ephys_sim';

% read in data simulated from python
sig0 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope0.00_%s.txt',fstem)),'ReadVariableNames',false);
sig01 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-0.10_%s.txt',fstem)),'ReadVariableNames',false);
sig02 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-0.20_%s.txt',fstem)),'ReadVariableNames',false);
sig03 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-0.30_%s.txt',fstem)),'ReadVariableNames',false);
sig04 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-0.40_%s.txt',fstem)),'ReadVariableNames',false);
sig05 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-0.50_%s.txt',fstem)),'ReadVariableNames',false);
sig06 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-0.60_%s.txt',fstem)),'ReadVariableNames',false);
sig07 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-0.70_%s.txt',fstem)),'ReadVariableNames',false);
sig08 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-0.80_%s.txt',fstem)),'ReadVariableNames',false);
sig09 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-0.90_%s.txt',fstem)),'ReadVariableNames',false);
sig1 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-1.00_%s.txt',fstem)),'ReadVariableNames',false);
sig11 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-1.10_%s.txt',fstem)),'ReadVariableNames',false);
sig12 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-1.20_%s.txt',fstem)),'ReadVariableNames',false);
sig13 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-1.30_%s.txt',fstem)),'ReadVariableNames',false);
sig14 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-1.40_%s.txt',fstem)),'ReadVariableNames',false);
sig15 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-1.50_%s.txt',fstem)),'ReadVariableNames',false);
sig16 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-1.60_%s.txt',fstem)),'ReadVariableNames',false);
sig17 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-1.70_%s.txt',fstem)),'ReadVariableNames',false);
sig18 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-1.80_%s.txt',fstem)),'ReadVariableNames',false);
sig19 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-1.90_%s.txt',fstem)),'ReadVariableNames',false);
sig2 = readtable(fullfile(datapath,sprintf('sim_neuralts_sig_slope-2.00_%s.txt',fstem)),'ReadVariableNames',false);

data4loop = {sig0.Var1,sig01.Var1,sig02.Var1,sig03.Var1,sig04.Var1, ...
    sig05.Var1,sig06.Var1,sig07.Var1,sig08.Var1,sig09.Var1, ...
    sig1.Var1,sig11.Var1,sig12.Var1,sig13.Var1,sig14.Var1,sig15.Var1, ...
    sig16.Var1,sig17.Var1,sig18.Var1,sig19.Var1,sig2.Var1};


for i = 1:length(data4loop)
    data2use = data4loop{i};
    data2use = detrend(data2use);

    % convolve with HRF
    data2use_BOLD = conv(data2use,hrf);

    % downsample
    data2use_BOLD_Down = downsample(data2use_BOLD,sampling_rate);
    data2use_BOLD_Down = data2use_BOLD_Down(1:fmri_tp);

    Heegraw(i,1) = bfn_mfin_ml(data2use,...
        'filter',Hfilter_type,'lb',lb_param,'ub',ub_param,'verbose',0);

    Hboldraw(i,1) = bfn_mfin_ml(data2use_BOLD_Down,...
        'filter',Hfilter_type,'lb',lb_param,'ub',ub_param,'verbose',0);

end % for i

legend_labels = {'1/f = 0','1/f = -0.1','1/f = -0.2','1/f = -0.3', ...
    '1/f = -0.4','1/f = -0.5','1/f = -0.6','1/f = -0.7','1/f = -0.8', ...
    '1/f = -0.9','1/f = -1','1/f = -1.1','1/f = -1.2','1/f = -1.3', ...
    '1/f = -1.4','1/f = -1.5','1/f = -1.6','1/f = -1.7','1/f = -1.8','1/f = -1.9','1/f = -2'};
legend_labels4table = {'0','-0.1','-0.2','-0.3', '-0.4','-0.5','-0.6', ...
    '-0.7','-0.8','-0.9','-1','-1.1','-1.2','-1.3','-1.4','-1.5','-1.6', ...
    '-1.7','-1.8','-1.9','-2'};

%%
if MAKE_PLOT
    disp('H from simulated BOLD when OOF is manipulated and gamma and alpha oscillations are added')
    figure; set(gcf,'color','w');
    subplot(4,1,1); scatter(1:length(legend_labels), Hbold_alphagamma'); grid on;
    hold on; lsline; % legend(legend_labels);
    set(gca,'XTick',1:length(legend_labels),'XTickLabels',legend_labels);
    rotateXLabels(gca,45); title('1/f + Gamma & Alpha (BOLD)');

    subplot(4,1,2); scatter(1:length(legend_labels), Hbold_alpha'); grid on;
    hold on; lsline; % legend(legend_labels);
    set(gca,'XTick',1:length(legend_labels),'XTickLabels',legend_labels);
    rotateXLabels(gca,45); title('1/f + Alpha (BOLD)');

    subplot(4,1,3); scatter(1:length(legend_labels), Hbold_gamma'); grid on;
    hold on; lsline; % legend(legend_labels);
    set(gca,'XTick',1:length(legend_labels),'XTickLabels',legend_labels);
    rotateXLabels(gca,45); title('1/f + Gamma (BOLD)');

    subplot(4,1,4); scatter(1:length(legend_labels), Heeg'); grid on;
    hold on; lsline; % legend(legend_labels);
    set(gca,'XTick',1:length(legend_labels),'XTickLabels',legend_labels);
    rotateXLabels(gca,45); title('1/f + Gamma & Alpha (EEG)');
end % if MAKE_PLOT

%%
% pack into results
result.oof = [0:0.1:2];
result.oof_legend_labels = legend_labels;
result.Hbold = Hboldraw';
result.Heeg = Heegraw';

%%
% write out to a file
varnames = {'legend_labels','oof','lfp','bold'};
tab2use = cell2table([[legend_labels4table'], ...
    num2cell([result.oof', result.Heeg', result.Hbold'])], ...
    'VariableNames',varnames);
writetable(tab2use,fullfile(datapath,sprintf('Hsim_%s.csv',fstem)),'FileType','text','delimiter',',');

end % function A_insilico_3_boldsim_neuralts_ooof3
