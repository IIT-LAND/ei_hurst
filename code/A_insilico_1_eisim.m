function result = A_insilico_1_eisim(EI_ratio, MAKE_PLOT)
% A_insilico_1_eisim - simulate LFP based on model from Gao, Peterson, & Voytek
%
% INPUT
%   EI_ratio = [2:0.2:6];
%   MAKE_PLOT = set to 1 to make plot
%
% Usage 
% EI_ratio = [2:0.2:6];
% MAKE_PLOT = 1;
% result = A_insilico_1_eisim(EI_ratio, MAKE_PLOT);
%

%%
outpath = '/Users/mlombardo/Dropbox/Manuscript/AIMS_Hurst/Sex/reproAnalysis/github_repo/data/gao_model';

%% set up simulation
dt = 0.001; %simulation time step
fs = 1/dt; %sampling rate
tk = 0:dt:1; %PSC kernel time vector
t = 0:dt:60*2; %simulation time vector

%spike train parameters
FR_E = 2;
FR_I = 5;
N_E = 8000;
N_I = 2000;

%ampa/gaba PSC kernels
Vr = -65;
Ee = 0;
Ei = -80;
AMPA_tau = [0.1 2]/1000;
GABA_tau = [0.5 10]/1000;
kA = syn_kernel(tk,AMPA_tau);
kG = syn_kernel(tk,GABA_tau);

%% simulate with different EI ratios
% simulate 128 regions
num_trs = 128;  
for i = 1:length(EI_ratio)
    
    leg_labels{i} = sprintf('1:%0.1f',EI_ratio(i));
    boost = EI_ratio(i)./((N_I*FR_I*sum(kG))/(N_E*FR_E*sum(kA)));
    disp(EI_ratio(i))   
    
    for tr = 1:num_trs        
        spk_E = pois_spikes(t(end)+tk(end), dt,N_E,FR_E);
        spk_I = pois_spikes(t(end)+tk(end), dt,N_I,FR_I);

        GE(:,i,tr) = conv(spk_E,kA, 'valid');
        GI(:,i,tr) = conv(spk_I,kG, 'valid')*boost;
        LFP_E(:,i,tr) = detrend(GE(:,i,tr),'constant')*(Ee-Vr);
        LFP_I(:,i,tr) = detrend(GI(:,i,tr),'constant')*(Ei-Vr);
        LFP(:,i,tr) = LFP_E(:,i,tr) + LFP_I(:,i,tr); 
    end % for tr
    
    % grab data to throw into function to compute H
    data2use = squeeze(LFP(:,i,:)); 
    % compute Hurst exponent 
    [H(:,i), NFC, FC] = bfn_mfin_ml(data2use,'filter','haar','lb',[-0.5,0],'ub',[1.5,10],'verbose',false);

end % for i


%% compute PSDs & fit
%oscillation mask
o_filt = ones(fs/2+1,1);
o_filt(9:13) = 1+4*gausswin(5); %alpha
o_filt(18:25) = 1+1.5*gausswin(8); %beta

o_mask = repmat(o_filt,[1, length(EI_ratio)]); %oscillation mask
for i = 1:size(LFP,3)
    %include oscillation
    PSD(:,:,i) = mPSD(LFP(:,:,i),fs,fs,fs/2,fs/2).*o_mask;
             
    rbf = robfitPSD(PSD(:,:,i),30:50,1); %robust fit    
    slope(:,i) = rbf(:,2);
end % for


if MAKE_PLOT
    % plot relationship between mean H across all regions and E:I ratio
    figure; set(gcf,'color','white');
    scatter(log2([1./EI_ratio])',fliplr(nanmean(H(:,:),1))');
    set(gca,'XTick',log2(fliplr(1./EI_ratio)),'XTickLabel',leg_labels);
    ylabel('Hurst exponent (H)'); xlabel('Excitation:Inhibition Ratio');
    rotateXLabels(gca,90);
    grid on;
    
    fname2save = 'EIratio_1to6.pdf';
    print(gcf,'-dpdf','-opengl','-r300',fname2save);
end % if MAKE_PLOT

%% write out result
for i = 1:num_trs
    varnames{i} = sprintf('channel_%03d',i);
end % for i
varnames = [{'EIlabels','EIratio'}, varnames];


tab2write = cell2table([leg_labels', num2cell([EI_ratio', slope])],'VariableNames',varnames);
fname2save = fullfile(outpath,'EIsim_oof.csv');
writetable(tab2write, fname2save, 'FileType','text','delimiter',',');

tab2write = cell2table([leg_labels', num2cell([EI_ratio', H'])],'VariableNames',varnames);
fname2save = fullfile(outpath,'EIsim_H.csv');
writetable(tab2write, fname2save, 'FileType','text','delimiter',',');
result = tab2write;

end % function A_insilico_1_eisim


%%
function kernel = syn_kernel(t, tau, type)
%function kernel = syn_kernel(t, tau, type)
% given a specific synaptic kernel type and time constant, this returns a
% time series of the kernel that spans the time defined (t) in seconds
%
% t: time vector in seconds (e.g. t=0:0.001:5)
% tau: t_decay or [t_rise t_decay] in seconds
% type: single- or double-exponential, or alpha synapse

if length(tau)==2
    type='2exp';
end % if length(tau)==2
switch type
    case 'exp'
        %single decay exponential (e^-t)
        kernel = exp(-t./tau(1));
        
    case '2exp'
        %double exponential -(e^t/t_rise) + e^t/t_decay
        if length(tau)==2
            % compute the normalization factor
            tpeak = tau(2)*tau(1)/(tau(2)-tau(1))*log(tau(2)/tau(1));
            normf = 1/(-exp(-tpeak./tau(1))+exp(-tpeak./tau(2)));            
            kernel = (-exp(-t./tau(1))+exp(-t./tau(2)))*normf;
        else
            kernel = [];
            disp('Need two time constants for double exponential.')
        end
                
    case 'alpha'
        %alpha kernel (t*e^t)
        kernel = (t/tau).*exp(-t./tau);        
end % switch type

end % function syn_kernel



%%
function discretized = pois_spikes(sim_t, dt, N_neu, FR)
%discretized = pois_spikes(sim_t, dt, N_neu, FR)
%   simulate population spiking of N neurons firing at FR each, return a
%   single spike train that is the total spiking

%mu parameter for exponential distribution
MU = 1./(N_neu*FR);

%draw ISI from exp RV 
ISI = exprnd(MU, [((sim_t+2)/MU) 1]); %pad 2s of sim time
spk_times = cumsum(ISI);
spk_times(spk_times>sim_t)=[];

%discretize
bins = (0:dt:sim_t)+dt/2; %make discretizing bins
discretized = hist(spk_times,bins);
end % function pois_spikes



%%
function [output, output_time, output_f] = stft(timestamp, data, Fs, winLen, stepLen, end_freq)
%function [output, output_time, output_f] = stft(timestamp, data, Fs, winLen, stepLen, end_freq)
% performs short-time windowed fourier transform on data, with a hamming
% window applied. Basically does the same thing as spectrogram.m
% 
% timestamp: leave blank []
% data: time series
% Fs: sampling frequency of the data
% winLen: window length in number of samples (use same number as Fs)
% stepLen: step length in number of samples (use about 1/50 as Fs (or 10 samples)
% end_freq (optional): cut of frequency of stft, default is fs/2

if isempty(timestamp)
   timestamp = linspace(1/Fs,size(data,1)/Fs, size(data,1)); 
end
%winLen=Fs;
max_iter = ceil((size(data,1)-winLen)/stepLen);
%preallocate output 
output_time = zeros(max_iter,1);
output_f = linspace(0,Fs-Fs/winLen, winLen);
if isempty(end_freq)
    end_ind = length(output_f);
else
    end_ind = find(output_f>end_freq,1)-1; %find index
end
output_f=output_f(1:end_ind);
    
output = zeros(end_ind, size(data,2), max_iter);
H = hamming(winLen);
HAM = repmat(H./sqrt(sum(H.^2)),1,size(data,2));

%%%stepping through
for i=1:max_iter
    %dc offset    
    cur_window = data((1:winLen)+(i-1)*stepLen,:)/sqrt(winLen);
    cur_window = (cur_window - repmat(mean(cur_window,1),winLen,1)).*HAM;         
    output_time(i)=timestamp(winLen+(i-1)*stepLen);    
    F = fft(cur_window);
    %use multitaper
    %F = pmtm(cur_window,4,winLen,Fs);
    output(:,:,i) = (F(1:end_ind,:));
    %output(:,:,i) = abs(F(1:end_ind,:));

    %output(:,i)=step_features(cur_window, FREQS);
end % for i
end % function stft



%%
function [PSD, Fa] = mPSD(data, Fs, winLen, stepLen, end_freq)
%function mPSD(data, Fs, winLen, stepLen, end_freq)
% computes median of stft, colloquially referred to as median welch.
% Basically it does the exact same thing as welch's method, but takes
% median instead of the mean, because median is robust to huge outliers and
% spectral data is non-uniformly distributed
%
% data: time series
% Fs: sampling frequency of the data
% winLen: window length in number of samples (use same number as Fs)
% stepLen: step length in number of samples (use about 1/50 as Fs (or 10 samples)
% end_freq (optional): cut of frequency of stft, default is fs/2
% PSD: returns median PSD
% Fa: returns frequency axis

if size(data,1)<size(data,2)
    data = data';
end

[F, Ft, Fa] = stft([], data, Fs, winLen, stepLen, end_freq);
PSD = median(abs(F).^2,3);
end % function mPSD


%%
function [fitparam, fiterr] = robfitPSD(PSD,freqs,dF)
%[fitparam, fiterr]= robfitPSD(PSD,freqs,dF)
% performs robust fit on PSD, over specifed frequency regions
%
%   PSD: can be 1,2 or 3D, must start from 0Hz
%       if 1D, will automatically rotate to have long axis as frequency
%       if 2D, needs to be []
%       if 3D, needs to be [X,Y,freq], where X & Y can be trial or time
%   freq: the frequency range to fit, should be at dF resolution
%   dF: is the frequency resolution of PSD [default 1]

if nargin == 2
    dF=1;
end
if length(size(PSD))==2
    %2D matrix, check if really 1D
    if any(size(PSD)==1)
       %really 1D, long axis is freq 
       if size(PSD,2)>size(PSD,1)
           PSD = PSD';
       end
    end
    PSD = permute(PSD,[2 3 1]); %expand into 3D array
end
[LX LY ~] = size(PSD);
fitparam = zeros(LX, LY, 2);
fiterr = zeros(LX, LY);
%fitInds = (round(freqs(1)/dF):round(freqs(end)/dF))+1;
fitInds = round(freqs/dF)+1; % find indices of PSD to fit over, corresponding to freqs
for xx = 1:LX
    for yy = 1:LY
        % perform fit and retrieve error
        [fitparam(xx,yy,:), temp] = robustfit(log10((fitInds-1)*dF),log10(squeeze(PSD(xx,yy,fitInds))));
        fiterr(xx,yy) = temp.ols_s;
    end
end
fitparam = squeeze(fitparam);
fiterr = squeeze(fiterr);
end % function robfitPSD


