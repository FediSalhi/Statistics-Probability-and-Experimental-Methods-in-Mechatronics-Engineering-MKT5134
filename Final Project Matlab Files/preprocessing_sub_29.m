%% specify your eeg data directory (EegMyDataDir) and temporary directory (TemDir)

EegMyDataDir = 'C:\Users\Designer\Desktop\Master''s Deg\Spring 2019 - 2020\Mekatronik Mühendisliðinde Ýstatistik Olasýlýk ve Deneysel Yöntemler\Project\Dataset\';
TemDir       = 'TempDir';
startup_bbci_toolbox('DataDir',EegMyDataDir,'TmpDir',TemDir); 
BTB.History = 0; % to aviod error for merging cnt

%% parameters
subdir_list = {'subject 26', 'subject 27', 'subject 28', 'subject 29'}; % subject
basename_list = {'motor_imagery1','mental_arithmetic1','motor_imagery2','mental_arithmetic1','motor_imagery3','mental_arithmetic1'}; % task type: motor imagery / recording session: 1 - 3
stimDef.eeg = {16,32; 'condition1','condition2'};

% load occular artifact-free eeg data (subject 29)
WorkingDir = 'C:\Users\Designer\Desktop\Master''s Deg\Spring 2019 - 2020\Mekatronik Mühendisliðinde Ýstatistik Olasýlýk ve Deneysel Yöntemler\Project\Matlab codes';
loadDir = fullfile(EegMyDataDir, subdir_list{4});
cd(loadDir);
load cnt; load mrk, load mnt; % load continous eeg signal (cnt), marker (mrk) and montage (mnt)
cd(WorkingDir)

%%
% merge cnts in each session
% for motor imagery: imag
% for mental arithmetic: ment
cnt_temp = cnt; mrk_temp = mrk; % save data temporarily
clear cnt mrk;
[cnt.imag, mrk.imag] = proc_appendCnt({cnt_temp{1}, cnt_temp{3}, cnt_temp{5}}, {mrk_temp{1}, mrk_temp{3}, mrk_temp{5}}); % merged motor imagery cnts
[cnt.ment, mrk.ment] = proc_appendCnt({cnt_temp{2}, cnt_temp{4}, cnt_temp{6}}, {mrk_temp{2}, mrk_temp{4}, mrk_temp{6}}); % merged mental arithmetic cnts



%% common average reference

% Average of all scalp channels (Common Average Reference)
% Re-referencing is achieved by creating an average of all scalp channels and subtracting the resulting 
% signal from each channel. After re-referencing, the overall electrical activity (amplitude) across all 
% channels will sum up to zero at each time point

cnt.imag = proc_commonAverageReference(cnt.imag);
cnt.ment = proc_commonAverageReference(cnt.ment);


%% separating left hand and right hand MI
left_hand_eeg   = zeros(1,32) ;
right_hand_eeg  = zeros(1,32) ;

for event_i=1:length(mrk.imag.time)-1
    chunk =  cnt.imag.x(round(mrk.imag.time(event_i)*0.2):round(mrk.imag.time(event_i+1)*0.2),:) ;
    if mrk.imag.event.desc(event_i) == 16
       right_hand_eeg = [right_hand_eeg; chunk];
    end
    
    if mrk.imag.event.desc(event_i) == 32
       left_hand_eeg = [left_hand_eeg; chunk];
    end
end

final_chunk = cnt.imag.x(round(mrk.imag.time(end)*0.2):end,:) ;
if mrk.imag.event.desc(event_i) == 16
       right_hand_eeg = [right_hand_eeg; chunk];
end
    
if mrk.imag.event.desc(event_i) == 32
       left_hand_eeg = [left_hand_eeg; chunk];
end


left_hand_eeg = left_hand_eeg(2:end,:);
right_hand_eeg = right_hand_eeg(2:end,:);

%% separating menatal arithmetic and resting state
mental_arithmetic_eeg   = zeros(1,32) ;
rest_eeg                = zeros(1,32) ;

for event_i=1:length(mrk.ment.time)-1
    chunk =  cnt.ment.x(round(mrk.ment.time(event_i)*0.2):round(mrk.ment.time(event_i+1)*0.2),:) ;
    if mrk.ment.event.desc(event_i) == 16
       rest_eeg = [rest_eeg; chunk];
    end
    
    if mrk.ment.event.desc(event_i) == 32
       mental_arithmetic_eeg = [mental_arithmetic_eeg; chunk];
    end
end

final_chunk = cnt.ment.x(round(mrk.ment.time(end)*0.2):end,:) ;
if mrk.ment.event.desc(event_i) == 16
       rest_eeg = [rest_eeg; chunk];
end
    
if mrk.ment.event.desc(event_i) == 32
       mental_arithmetic_eeg = [mental_arithmetic_eeg; chunk];
end


mental_arithmetic_eeg = mental_arithmetic_eeg(2:end,:);
rest_eeg = rest_eeg(2:end,:);

%% outliers detection and removal
%plot distribution
figure(1), clf, hold on
subplot(2,2,1)
histogram(right_hand_eeg(:,1)),
title('Right Hand Motor Imagery'),
xlabel('Channel 1 - Voltage (mV)')
ylabel('Freqeuncy');

subplot(2,2,2)
histogram(left_hand_eeg(:,1)),
title('Left Hand Motor Imagery'),
xlabel('Channel 1 - Voltage (mV)')
ylabel('Freqeuncy');

subplot(2,2,3)
histogram(mental_arithmetic_eeg(:,1)),
title('Mental Arithmetic '),
xlabel('Channel 1 - Voltage (mV)')
ylabel('Freqeuncy');

subplot(2,2,4)
histogram(rest_eeg(:,1)),
title('Rest State'),
xlabel('Channel 1 - Voltage (mV)')
ylabel('Freqeuncy');

left_hand_eeg_no_outliers         = rmoutliers(left_hand_eeg,'quartiles');
right_hand_eeg_no_outliers        = rmoutliers(right_hand_eeg,'quartiles');
mental_arithmetic_eeg_no_outliers = rmoutliers(mental_arithmetic_eeg,'quartiles');
rest_eeg_no_outliers              = rmoutliers(rest_eeg,'quartiles');



%% plot signals before and after outliers removal
figure(2),
subplot(1,2,1),
boxplot([right_hand_eeg(1:165000,1), left_hand_eeg(1:165000,1), mental_arithmetic_eeg(1:165000,1), rest_eeg(1:165000,1)], 'Labels', {'Right hand MI', 'Left hand MI', 'MA', 'Rest state'})
title('Before removing outliers')
xlabel('Event')
%ylim([-350 250
ylabel('Channel-1 volatage (mV)')

subplot(1,2,2),
boxplot([right_hand_eeg_no_outliers(1:107648,1), left_hand_eeg_no_outliers(1:107648,1), mental_arithmetic_eeg_no_outliers(1:107648,1), rest_eeg_no_outliers(1:107648,1)], 'Labels', {'Right hand MI', 'Left hand MI', 'MA', 'Rest state'})
title('After removing outliers')
xlabel('Event')
ylabel('Channel-1 volatage (mV)')

%% bandpass filter
% In order to reduce the effect of the voltage from the galvanic skin response across the head, 
% the EEG may use a filter to attenuate frequencies below 5 Hz.
% A low pass filter is needed to reject noise at higher frequencies than the EEG signals.

fs = 200; %Hz
fpass = [5 50] ; %hz
right_hand_eeg_no_outliers_filtered         = bandpass(right_hand_eeg_no_outliers,fpass,fs);
left_hand_eeg_no_outliers_filtered          = bandpass(left_hand_eeg_no_outliers,fpass,fs);
mental_arithmetic_eeg_no_outliers_filtered  = bandpass(mental_arithmetic_eeg_no_outliers,fpass,fs);
rest_eeg_no_outliers_filtered               = bandpass(rest_eeg_no_outliers,fpass,fs);

%% Principle component analysis - PCA
[coeff_right_hand_MI,score,latent,tsquared,explained_right_hand_MI,mu] = pca(right_hand_eeg_no_outliers_filtered);
[coeff_left_hand_MI,score,latent,tsquared,explained_left_hand_MI,mu]   = pca(left_hand_eeg_no_outliers_filtered);
[coeff_MA,score,latent,tsquared,explained_MA,mu]                       = pca(mental_arithmetic_eeg_no_outliers_filtered);
[coeff_rest,score,latent,tsquared,explained_rest,mu]                   = pca(rest_eeg_no_outliers_filtered);

%% linear combinations
%right hand MI
[r1, c1] = size(right_hand_eeg_no_outliers_filtered);
right_hand_eeg_no_outliers_filtered_combined = zeros(r1,1);

for i=1:32
    component=coeff_right_hand_MI(i)*right_hand_eeg_no_outliers_filtered(:,i);
    right_hand_eeg_no_outliers_filtered_combined = right_hand_eeg_no_outliers_filtered_combined + component ;
end

%left hand MI
[r2, c2] = size(left_hand_eeg_no_outliers_filtered);
left_hand_eeg_no_outliers_filtered_combined = zeros(r2,1);

for i=1:32
    component=coeff_left_hand_MI(i)*left_hand_eeg_no_outliers_filtered(:,i);
    left_hand_eeg_no_outliers_filtered_combined = left_hand_eeg_no_outliers_filtered_combined + component ;
end

%MA
[r3, c3] = size(mental_arithmetic_eeg_no_outliers_filtered);
mental_arithmetic_eeg_no_outliers_filtered_combined = zeros(r3,1);

for i=1:32
    component=coeff_MA(i)*mental_arithmetic_eeg_no_outliers_filtered(:,i);
    mental_arithmetic_eeg_no_outliers_filtered_combined = mental_arithmetic_eeg_no_outliers_filtered_combined + component ;
end

%rest
[r4, c4] = size(rest_eeg_no_outliers_filtered);
rest_eeg_no_outliers_filtered_combined = zeros(r4,1);

for i=1:32
    component=coeff_rest(i)*rest_eeg_no_outliers_filtered(:,i);
    rest_eeg_no_outliers_filtered_combined = rest_eeg_no_outliers_filtered_combined + component ;
end

%% computing average power
% power spectrum density
right_hand_eeg_no_outliers_filtered_combined_Hpsd = dspdata.psd(abs(right_hand_eeg_no_outliers_filtered_combined),'Fs',fs); %power spectrum density estimation 
left_hand_eeg_no_outliers_filtered_combined_Hpsd = dspdata.psd(abs(left_hand_eeg_no_outliers_filtered_combined),'Fs',fs);
mental_arithmetic_eeg_no_outliers_filtered_combined_Hpsd = dspdata.psd(abs(mental_arithmetic_eeg_no_outliers_filtered_combined),'Fs',fs);
rest_eeg_no_outliers_filtered_combined_Hpsd = dspdata.psd(abs(rest_eeg_no_outliers_filtered_combined),'Fs',fs);

% average power
avrRH = avgpower(right_hand_eeg_no_outliers_filtered_combined_Hpsd)
avrLH = avgpower(left_hand_eeg_no_outliers_filtered_combined_Hpsd)
avrMA = avgpower(mental_arithmetic_eeg_no_outliers_filtered_combined_Hpsd)
avrRes= avgpower(rest_eeg_no_outliers_filtered_combined_Hpsd)
