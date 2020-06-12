%% C_3a_AIMS_Hurst_PLS.m

datapath = '/Users/mlombardo/Dropbox/Manuscripts/AIMS_Hurst_Sex/reproAnalysis/data/mrc_aims/Hurst_dwt';
phenopath = '/Users/mlombardo/Dropbox/Manuscripts/AIMS_Hurst_Sex/reproAnalysis/pheno';
plsrespath = '/Users/mlombardo/Dropbox/Manuscripts/AIMS_Hurst_Sex/reproAnalysis/pls_results';

pdata = readtable(fullfile(phenopath,'tidy_data.csv'));
mask = pdata.use_subs==1;
sublist = pdata.sub_id(mask);

tmp_data = pdata(mask,:);

td_mask = ismember(tmp_data.Diagnosis,'TD');
asd_mask = ~td_mask;
male_mask = ismember(tmp_data.Sex,'M');
female_mask = ~male_mask;

ngrp = 4;
design_mat = zeros(size(tmp_data,1),ngrp);
design_mat(td_mask & male_mask,1) = 1;
design_mat(asd_mask & male_mask,2) = 1;
design_mat(td_mask & female_mask,3) = 1;
design_mat(asd_mask & female_mask,4) = 1;

contrast_mat = zeros(size(design_mat,1),3);
contrast_mat(:,1) = design_mat * [1 -1 1 -1]';
contrast_mat(:,2) = design_mat * [1 1 -1 -1]';
contrast_mat(:,3) = design_mat * [1 -1 -1 1]';

nreg = 180;

hdata = zeros(length(sublist),nreg);

for isub = 1:length(sublist)
    fname = fullfile(datapath,sprintf('%d__Hurst.csv',sublist(isub)));
    disp(fname);
    tmpdata = readtable(fname,'ReadVariableNames',false);
    hdata(isub,:) = table2array(tmpdata);
end % for i

X = {hdata};
Y = contrast_mat;


%% Analysis parameters

%Number of boot and perm
nperm = 10000;

% Number of subjects
num_subj = [];
for N=1:numel(X)
    num_subj(N) = size(X{N},1);
end

num_cond = 1; 

% specify option structure
option.method = 3; %[1] | 2 | 3 | 4 | 5 | 6
option.num_perm = nperm; %( single non-negative integer )
option.is_struct = 0;%[0] | 1
option.num_split = 0; %( single non-negative integer )
option.num_boot = nperm; %( single non-negative integer )
option.clim = 95; %( [95] single number between 0 and 100 )
option.stacked_behavdata = Y;
option.cormode = 0; %[0] | 2 | 4 | 6
option.boot_type = 'strat'; %['strat'] | 'nonstrat'


%% Run PLS on ALL CONTRASTS
result = pls_analysis(X, num_subj, num_cond, option);

% Compute percentage of cross-block covariance
result.crossblockCovPercent = result.s.^2/sum(result.s.^2);

% fix p-values
result.perm_result.sprob = (result.perm_result.sp+1)./(result.perm_result.num_perm+1);

% find the significant LVs
result.sigLVs = find(result.perm_result.sprob<=0.05);
result.sigLVs_pvals = result.perm_result.sprob(result.sigLVs);
disp(sprintf('Significant LVs: LV %d',result.sigLVs));

% compute brain BSR
for i = 1:size(result.boot_result.compare_u,2)
    result.brain_bsr(:,i) = result.boot_result.compare_u(:,i)./result.boot_result.u_se(:,i);
end

% make parcels for this df
GlasserRegions = 180;
Parcels = cell(GlasserRegions, 1);
for reg=1:GlasserRegions
    region = cellstr(sprintf('parcel_%03s', string(reg)));
    Parcels(reg) = region;
end

result.brainreg_names = Parcels';

% save results to file
save(fullfile(plsrespath,'pls_ALL_H.mat'), 'result');

% output brain BSR
for i = 1:length(result.sigLVs)
    bsr_names{i} = sprintf('BSR_LV%d',result.sigLVs(i));
end
tab2write = cell2table([result.brainreg_names', num2cell([[1:nreg]', result.brain_bsr(:,result.sigLVs)])],'VariableNames',[{'brainreg','parcel_index'},bsr_names]);
writetable(tab2write, fullfile(plsrespath,'pls_ALL_H_brainBSR4plotting.csv'),'FileType','text','Delimiter',',');

% compute bootstrap CIs
cis2use = {[2.5,97.5]};
for i = 1:length(result.sigLVs)
    LVnum = result.sigLVs(i);
    for j = 1:length(cis2use)
        ci_bounds = cis2use{j};

        dx_idx = 1;
        sex_idx = 2;
        dxxsex_idx = 3;

        dx_bootres = squeeze(result.boot_result.distrib(dx_idx,LVnum,:));
        sex_bootres = squeeze(result.boot_result.distrib(sex_idx,LVnum,:));
        dxxsex_bootres = squeeze(result.boot_result.distrib(dxxsex_idx,LVnum,:));

        dx_ci = prctile(dx_bootres',ci_bounds)';
        sex_ci = prctile(sex_bootres',ci_bounds)';
        dxxsex_ci = prctile(dxxsex_bootres',ci_bounds)';

        dx_corr = result.boot_result.orig_corr(dx_idx,LVnum);
        sex_corr = result.boot_result.orig_corr(sex_idx,LVnum);
        dxxsex_corr = result.boot_result.orig_corr(dxxsex_idx,LVnum);

        dx_data = [dx_corr dx_ci'];
        sex_data = [sex_corr sex_ci'];
        dxxsex_data = [dxxsex_corr dxxsex_ci'];

        dx_labels = cellstr(repmat({'Dx'},1,1));
        sex_labels = cellstr(repmat({'Sex'},1,1));
        dxxsex_labels = cellstr(repmat({'Dx_x_Sex'},1,1));
        effect_labels = [dx_labels';sex_labels';dxxsex_labels'];

        results2use = [dx_data;sex_data;dxxsex_data];

        tab2write = cell2table([effect_labels, num2cell(results2use)], ...
            'VariableNames',{'Effect','corr','lo_lim','up_lim'});

        file2save = fullfile(plsrespath,sprintf('pls_ALL_H_bootCI4plotting_LV%d_ci%d.csv',LVnum,ci_bounds(2)-ci_bounds(1)));
        writetable(tab2write,file2save,'FileType','text','delimiter',',');
    end % for j
end % for i


%% Run PLS on SEXbyDX CONTRAST
clear result;
option.stacked_behavdata = Y(:,3);
result = pls_analysis(X, num_subj, num_cond, option);

% Compute percentage of cross-block covariance
result.crossblockCovPercent = result.s.^2/sum(result.s.^2);

% fix p-values
result.perm_result.sprob = (result.perm_result.sp+1)./(result.perm_result.num_perm+1);

% find the significant LVs
result.sigLVs = find(result.perm_result.sprob<=0.05);
result.sigLVs_pvals = result.perm_result.sprob(result.sigLVs);
disp(sprintf('Significant LVs: LV %d',result.sigLVs));

% compute brain BSR
for i = 1:size(result.boot_result.compare_u,2)
    result.brain_bsr(:,i) = result.boot_result.compare_u(:,i)./result.boot_result.u_se(:,i);
end

% make parcels for this df
GlasserRegions = 180;
Parcels = cell(GlasserRegions, 1);
for reg=1:GlasserRegions
    region = cellstr(sprintf('parcel_%03s', string(reg)));
    Parcels(reg) = region;
end

result.brainreg_names = Parcels';

% save results to file
save(fullfile(plsrespath,'pls_SEXbyDX_H.mat'), 'result');

% output brain BSR
for i = 1:length(result.sigLVs)
    bsr_names{i} = sprintf('BSR_LV%d',result.sigLVs(i));
end
tab2write = cell2table([result.brainreg_names', num2cell([[1:nreg]', result.brain_bsr(:,result.sigLVs)])],'VariableNames',[{'brainreg','parcel_index'},bsr_names]);
writetable(tab2write, fullfile(plsrespath,'pls_SEXbyDX_H_brainBSR4plotting.csv'),'FileType','text','Delimiter',',');

% compute bootstrap CIs
cis2use = {[2.5,97.5]};
for i = 1:length(result.sigLVs)
    LVnum = result.sigLVs(i);
    for j = 1:length(cis2use)
        ci_bounds = cis2use{j};

        dxxsex_bootres = squeeze(result.boot_result.distrib(:,LVnum,:));
        dxxsex_ci = prctile(dxxsex_bootres',ci_bounds)';
        dxxsex_corr = result.boot_result.orig_corr(:,LVnum);

        dxxsex_data = [dxxsex_corr dxxsex_ci'];
        dxxsex_labels = cellstr(repmat({'Dx_x_Sex'},1,1));
        effect_labels = [dxxsex_labels'];

        results2use = [dxxsex_data];

        tab2write = cell2table([effect_labels, num2cell(results2use)], ...
            'VariableNames',{'Effect','corr','lo_lim','up_lim'});

        file2save = fullfile(plsrespath,sprintf('pls_SEXbyDX_H_bootCI4plotting_LV%d_ci%d.csv',LVnum,ci_bounds(2)-ci_bounds(1)));
        writetable(tab2write,file2save,'FileType','text','delimiter',',');

        boot_labels = cell(nperm,1);
        for iperm = 1:nperm
            boot_labels{iperm} = sprintf('boot_%03d',iperm);
        end
        tab2write = cell2table([[{'actual_corr'};boot_labels], num2cell(dxxsex_bootres)],'VariableNames',{'boot_labels','boot_res'});
        file2save = fullfile(plsrespath,sprintf('pls_SEXbyDX_H_bootres4plotting_LV%d_ci%d.csv',LVnum,ci_bounds(2)-ci_bounds(1)));
        writetable(tab2write,file2save,'FileType','text','delimiter',',');
    end % for j
end % for i



%% Run PLS on SEX CONTRAST
clear result;
option.stacked_behavdata = Y(:,2);
result = pls_analysis(X, num_subj, num_cond, option);

% Compute percentage of cross-block covariance
result.crossblockCovPercent = result.s.^2/sum(result.s.^2);

% fix p-values
result.perm_result.sprob = (result.perm_result.sp+1)./(result.perm_result.num_perm+1);

% find the significant LVs
result.sigLVs = find(result.perm_result.sprob<=0.05);
result.sigLVs_pvals = result.perm_result.sprob(result.sigLVs);
disp(sprintf('Significant LVs: LV %d',result.sigLVs));

% compute brain BSR
for i = 1:size(result.boot_result.compare_u,2)
    result.brain_bsr(:,i) = result.boot_result.compare_u(:,i)./result.boot_result.u_se(:,i);
end

% make parcels for this df
GlasserRegions = 180;
Parcels = cell(GlasserRegions, 1);
for reg = 1:GlasserRegions
    region = cellstr(sprintf('parcel_%03s', string(reg)));
    Parcels(reg) = region;
end % for reg

result.brainreg_names = Parcels';

% save results to file
save(fullfile(plsrespath,'pls_SEX_H.mat'), 'result');

% output brain BSR
for i = 1:length(result.sigLVs)
    bsr_names{i} = sprintf('BSR_LV%d',result.sigLVs(i));
end
tab2write = cell2table([result.brainreg_names', num2cell([[1:nreg]', result.brain_bsr(:,result.sigLVs)])],'VariableNames',[{'brainreg','parcel_index'},bsr_names]);
writetable(tab2write, fullfile(plsrespath,'pls_SEX_H_brainBSR4plotting.csv'),'FileType','text','Delimiter',',');

% compute bootstrap CIs
cis2use = {[2.5,97.5]};
for i = 1:length(result.sigLVs)
    LVnum = result.sigLVs(i);
    for j = 1:length(cis2use)
        ci_bounds = cis2use{j};
        
        sex_bootres = squeeze(result.boot_result.distrib(:,LVnum,:));
        sex_ci = prctile(sex_bootres',ci_bounds)';
        sex_corr = result.boot_result.orig_corr(:,LVnum);
        sex_data = [sex_corr sex_ci'];
        sex_labels = cellstr(repmat({'Sex'},1,1));
        effect_labels = [sex_labels'];

        results2use = [sex_data];

        tab2write = cell2table([effect_labels, num2cell(results2use)], ...
            'VariableNames',{'Effect','corr','lo_lim','up_lim'});

        file2save = fullfile(plsrespath,sprintf('pls_SEX_H_bootCI4plotting_LV%d_ci%d.csv',LVnum,ci_bounds(2)-ci_bounds(1)));
        writetable(tab2write,file2save,'FileType','text','delimiter',',');

        boot_labels = cell(nperm,1);
        for iperm = 1:nperm
            boot_labels{iperm} = sprintf('boot_%03d',iperm);
        end % for iperm
        tab2write = cell2table([[{'actual_corr'};boot_labels], num2cell(sex_bootres)],'VariableNames',{'boot_labels','boot_res'});
        file2save = fullfile(plsrespath,sprintf('pls_SEX_H_bootres4plotting_LV%d_ci%d.csv',LVnum,ci_bounds(2)-ci_bounds(1)));
        writetable(tab2write,file2save,'FileType','text','delimiter',',');
    end % for j
end % for i



%% Run PLS on DX CONTRAST
clear result;
option.stacked_behavdata = Y(:,1);
result = pls_analysis(X, num_subj, num_cond, option);

% Compute percentage of cross-block covariance
result.crossblockCovPercent = result.s.^2/sum(result.s.^2);

% fix p-values
result.perm_result.sprob = (result.perm_result.sp+1)./(result.perm_result.num_perm+1);

% find the significant LVs
result.sigLVs = find(result.perm_result.sprob<=0.05);
if isempty(result.sigLVs)
    result.sigLVs_pvals = NaN;
    disp('No significant LVs');
else
    result.sigLVs_pvals = result.perm_result.sprob(result.sigLVs);
    disp(sprintf('Significant LVs: LV %d',result.sigLVs));
end

% compute brain BSR
for i = 1:size(result.boot_result.compare_u,2)
    result.brain_bsr(:,i) = result.boot_result.compare_u(:,i)./result.boot_result.u_se(:,i);
end

% make parcels for this df
GlasserRegions = 180;
Parcels = cell(GlasserRegions, 1);
for reg=1:GlasserRegions
    region = cellstr(sprintf('parcel_%03s', string(reg)));
    Parcels(reg) = region;
end

result.brainreg_names = Parcels';

% save results to file
save(fullfile(plsrespath,'pls_DX_H.mat'), 'result');

% output brain BSR
if ~isempty(result.sigLVs)
    for i = 1:length(result.sigLVs)
        bsr_names{i} = sprintf('BSR_LV%d',result.sigLVs(i));
    end
    tab2write = cell2table([result.brainreg_names', num2cell([[1:nreg]', result.brain_bsr(:,result.sigLVs)])],'VariableNames',[{'brainreg','parcel_index'},bsr_names]);
    writetable(tab2write, fullfile(plsrespath,'pls_DX_H_brainBSR4plotting.csv'),'FileType','text','Delimiter',',');


    % compute bootstrap CIs
    cis2use = {[2.5,97.5]};
    for i = 1:length(result.sigLVs)
        LVnum = result.sigLVs(i);
        for j = 1:length(cis2use)
            ci_bounds = cis2use{j};

            dx_bootres = squeeze(result.boot_result.distrib(:,LVnum,:));
            dx_ci = prctile(dx_bootres',ci_bounds)';
            dx_corr = result.boot_result.orig_corr(:,LVnum);
            dx_data = [dx_corr dx_ci'];
            dx_labels = cellstr(repmat({'Dx'},1,1));
            effect_labels = [dx_labels'];

            results2use = [dx_data];

            tab2write = cell2table([effect_labels, num2cell(results2use)], ...
                'VariableNames',{'Effect','corr','lo_lim','up_lim'});

            file2save = fullfile(plsrespath,sprintf('pls_DX_H_bootCI4plotting_LV%d_ci%d.csv',LVnum,ci_bounds(2)-ci_bounds(1)));
            writetable(tab2write,file2save,'FileType','text','delimiter',',');

            boot_labels = cell(nperm,1);
            for iperm = 1:nperm
                boot_labels{iperm} = sprintf('boot_%03d',iperm);
            end % for iperm
            tab2write = cell2table([[{'actual_corr'};boot_labels], num2cell(dx_bootres)],'VariableNames',{'boot_labels','boot_res'});
            file2save = fullfile(plsrespath,sprintf('pls_DX_H_bootres4plotting_LV%d_ci%d.csv',LVnum,ci_bounds(2)-ci_bounds(1)));
            writetable(tab2write,file2save,'FileType','text','delimiter',',');
        end % for j
    end % for i

end % if ~isempty


