%Combining parcellate.m and bfn_mfin_ml.m in MRC subjets
%Pipeline runs on resampled Glasser Temples to native space
%Assuming 3dresample with Erest_pp as -master was run and result was saved in a
%FINALHURSTS dir,
%and its the only one in the dir


%%

%Paths
rootpath = '/home/stavros/Documents/ABIDE_1/MRC_AIMS/rsfMRI/raw_data/TEST_compilation_script';


%Parellate args
MEANCENTER = 1;
nreg = 180;

%Non-fractal-master/m/bfn_mfin_ml.m args
lb = [-0.5,0];
ub = [1.5,10];

%Get subdir list - every subjets has a 0 in their id
d = dir(fullfile(rootpath, '*0*'));

%Loop

for i=1:length(d)
    subname = d(i).name;
    disp(subname);
    subpath = fullfile(rootpath, subname);
    disp(subpath);

    %csv path
    csvpath = fullfile(subpath, 'preproc_2','CSV');
    mkdir(csvpath);

    %parcellate args
    atlasfile = fullfile(subpath,'preproc_2', 'FINALHURSTS', 'Erest_TEMP.nii.gz');
    datafile = fullfile(subpath, 'preproc_2', 'Erest_pp.nii.gz');
    % saving time-series
    fname2save = fullfile(csvpath, strcat(subname,'_MMP_in_MNI_symmetrical_1.csv'));


    %calls parcellate.m
    result = C_1b_parcellate(atlasfile,datafile,fname2save,MEANCENTER, nreg);


   %import csv make it an array
    table = readtable(fname2save);
    disp(size(table));
    ArrayTable = table2array(table);

    % Identify if there are any regions that the atlas is not covering
    colswnan = find(sum(isnan(ArrayTable),1)~=0);


    %feed array into bfn_mfin_ml.m
    [H, nfcor, fcor] = bfn_mfin_ml(ArrayTable, 'filter', 'Haar', 'lb', lb, 'ub', ub);

    %if there are any columns with NAN replace them in ArrayTable with NAN
	if isempty(colswnan)
		H(colswnan) = NaN;
        nfcor(colswnan) = NaN;
        fcor(colswnan) = NaN;
    end



    %write everything into csv for each subject
    H_name = fullfile(csvpath, strcat(subname, '_H.csv'));
    fcor_name = fullfile(csvpath, strcat(subname, '_fcor.csv'));
    nfcor_name = fullfile(csvpath, strcat(subname, '_nfcor.csv'));


    csvwrite(H_name, H);
    csvwrite(fcor_name, fcor);
    csvwrite(nfcor_name, nfcor);


    clear H fcor nfcor datafile atlasfile fname2save H_name fcor_name nfcor_name att

end % ends for
