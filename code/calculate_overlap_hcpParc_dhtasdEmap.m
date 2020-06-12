% calculate_overlap_hcpParc_dhtasdEmap.m

addpath /Users/mlombardo/Dropbox/matlab/spm12

rootpath = '/Users/mlombardo/Dropbox/Manuscripts/AIMS_Hurst_Sex/reproAnalysis/dht_DE_asdEgenes';
codepath = '/Users/mlombardo/Dropbox/Manuscripts/AIMS_Hurst_Sex/reproAnalysis/code';

wb_path = fullfile(rootpath,'whole_brain_Satterstrom_Velmeshev_ExGenes_DHT');

mask = spm_read_vols(spm_vol(fullfile(wb_path,'mask.nii'))); mask = mask>0;

fdr_thresh = 3.148288;
tmap = spm_read_vols(spm_vol(fullfile(wb_path,'spmT_0001.nii')));
fdr_map = tmap>=fdr_thresh;

hcp_annot = readtable(fullfile(codepath,'GlasserHCP_annot.txt'));
hcp_annot.percent_overlap = zeros(size(hcp_annot,1),1);
hcp = spm_read_vols(spm_vol(fullfile(wb_path,'MMP_in_MNI_symmetrical_1_resamp.nii')));
nregs = 180;
hcp_regs = 1:nregs;

for i = 1:nregs
    reg2use = mask & ismember(hcp,i); 
    ntotal = sum(reg2use(:));
    overlap_map = reg2use & fdr_map;
    noverlap = sum(overlap_map(:));
    hcp_annot.percent_overlap(i) = (noverlap/ntotal)*100;
end % for i

fname2save = fullfile(wb_path,'hcp_fdr01_percent_overlap.csv');
writetable(hcp_annot,fname2save,'FileType','text','delimiter',',');