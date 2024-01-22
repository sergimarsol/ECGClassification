% Biomedical Engineering: Biomedical Signals and Instrumentation
% Sergi Marsol Torrent
% December 2022

%% Import data & plot one from each type (bandpower and periodogram)
root_cd = cd("ECG_normal/");
dir_info = struct2cell(dir())';
dir_true = cell2mat(dir_info(:,5));
file_true = ~dir_true;
files = dir_info(file_true,1);
load(files{2});
AN=val;
load(files{8});
AA=val;
load(files{16});
AO=val;
load(files{44});
AW=val;

%plots for each
figure(1), plot(AA)
figure(2), plot(AN)
figure(3), plot(AO)
figure(4), plot(AW)

cd('..')
cd('ECG_noisy/') 
root_cd = cd();
dir_info = struct2cell(dir())';
dir_true = cell2mat(dir_info(:,5));
file_true = ~dir_true;
files2 = dir_info(file_true,1);

load(files2{2});
BN=newval;
load(files2{5});
BN2=newval;
load(files2{10});
BN3=newval;
load(files2{15});
BN4=newval;
figure(5), plot(BN)

cd('..')
%periodogram for each
fs=300
figure(6), periodogram(AN,[],[],fs);
figure(7), periodogram(AA,[],[],fs);
figure(8), periodogram(AO,[],[],fs);
figure(9), periodogram(AW,[],[],fs);
figure(10), periodogram(BN,[],[],fs);
figure(11), periodogram(BN2,[],[],fs);
figure(12), periodogram(BN3,[],[],fs);
figure(13), periodogram(BN4,[],[],fs);

%bandpower for each
freqrange=[49 51];
b1=bandpower(AN,fs,freqrange)
b2=bandpower(AA,fs,freqrange)
b3=bandpower(AO,fs,freqrange)
b4=bandpower(AW,fs,freqrange)
b5=bandpower(BN,fs,freqrange)


%% Compute bandpower for all A and B
root_cd = cd("training2017");
dir_info = struct2cell(dir())';
dir_true = cell2mat(dir_info(:,5));
file_true = ~dir_true;
files = dir_info(file_true,1);
vectorA=zeros((numel(files))/2,1);
fs=300;
freqrange=[49 51];
A = zeros((numel(files)/2),18286);
vectorA_re = vectorA;
for i=1:(numel(files)/2)
    load(files{i*2})
    A(i,1:numel(val))=val;
    vectorA(i)=bandpower(val,fs,freqrange);
end

cd('..')
cd('ECG_noisy') 
root_cd = cd();
dir_info = struct2cell(dir())';
dir_true = cell2mat(dir_info(:,5));
file_true = ~dir_true;
files2 = dir_info(file_true,1);
vectorB=zeros(numel(files2),1);
vectorB_re = vectorB;
B = zeros(numel(files2),18286);
for i=1:numel(files2)
    load(files2{i})
    B(i,1:numel(newval))=newval;
    vectorB(i)=bandpower(newval,fs,freqrange);
end
cd('..')

%% Compute bandpower after resampling and normalization
A_re = zeros(8528,12191);
B_re = zeros(4999,12191);

for i=1:8528
    A_re(i,:)=resample(A(i,:),2,3);
    A_re(i,:)=A_re(i,:)./max(A_re(i,:));
    vectorA_re(i) =bandpower(A_re(i,:),200, freqrange);
end

for i=1:4999
    B_re(i,:)=resample(B(i,:),2,3);
    B_re(i,:)=B_re(i,:)./max(B_re(i,:));
    vectorB_re(i) = bandpower(B_re(i,:),200,freqrange);
end

vectorA_re = (vectorA_re - mean(vectorA_re))/std(vectorA_re);
vectorB_re = (vectorB_re - mean(vectorB_re))/std(vectorB_re);

%% Pan-Tomkins algorithm to detect QRS complexes in the EKG signal
%Step 1: preprocessing of the signal with filters, derivatives and
%integrations

filtA = zeros(8528,12191);
filtA2 = zeros(8528,12191);
derA = zeros(8528,12191);
derA2 = zeros(8528,12191);
derA3 = zeros(8528,12191);
sqA = zeros(8528,12191);
intA = zeros(8528,12191);

a0=[1 -2 1];
b0=[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a00=[1 -1];
b00=[-1/32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1/32];
a=1;
b1=[1 0 -1];
b2=[1 0 -2 0 1];
b3=[0.125 0.125 0.125 0.125 0.125 0.125 0.125 0.125];

for i=1:8528
    filtA(i,:) = filter(b0,a0,A_re(i,:));
    filtA(i,:) = filtA(i,:)./max(filtA(i,:));
    filtA2(i,:) = filter(b00,a00,filtA(i,:));
    filtA2(i,:) = filtA2(i,:)./max(filtA2(i,:));
    derA(i,:) = filter(b1,a,filtA2(i,:));
    derA(i,:) = derA(i,:)./max(derA(i,:));
    derA2(i,:) = filter(b2,a,filtA2(i,:));
    derA2(i,:) = derA2(i,:)./max(derA2(i,:));
    derA3(i,:) = 1.3*derA(i,:) + 1.1*derA2(i,:);
    derA3(i,:) = derA3(i,:)./max(derA3(i,:));
    sqA(i,:) = derA3(i,:).^2;
    sqA(i,:) = sqA(i,:)./max(sqA(i,:));
    intA(i,:) = filter(b3,a,sqA(i,:));
    intA(i,:) = intA(i,:)./max(intA(i,:));
end

%% plot of the previous steps

figure(14)
subplot(6,1,1)
plot(A_re(1,:))
xlim([4000 5000])
title("Original signal")
xlabel("Samples (n)")
ylabel("Amplitude")
subplot(6,1,2)
plot(filtA2(1,:))
xlim([4000 5000])
title("After PassBand Filter")
xlabel("Samples (n)")
ylabel("Amplitude")
subplot(6,1,3)
plot(derA(1,:))
xlim([4000 5000])
title("First derivative filter")
xlabel("Samples (n)")
ylabel("Amplitude")
subplot(6,1,4)
plot(derA2(1,:))
xlim([4000 5000])
title("Second derivative filter")
xlabel("Samples (n)")
ylabel("Amplitude")
subplot(6,1,5)
plot(derA3(1,:))
xlim([4000 5000])
title("Sum of derivatives")
xlabel("Samples (n)")
ylabel("Amplitude")
subplot(6,1,6)
plot(intA(1,:))
xlim([4000 5000])
title("Square integrated signal")
xlabel("Samples (n)")
ylabel("Amplitude")


%% Step 2: local maxima
% find the local maxima quantity
num_pks = zeros(1,8528);

for i=1:8528
    [pks, locs] = findpeaks(intA(i,:), "MinPeakDistance", 40);
    num_pks(i) = numel(pks);
end

max_pks =max(num_pks);

%% Find the local maxima of each (fiducial mark)
pks = zeros(max_pks,8528);
locs = zeros(max_pks,8528);

for i=1:8528
    [pks(1:num_pks(i),i),locs(1:num_pks(i),i)]=findpeaks(intA(i,:), "MinPeakDistance", 40);
end
% the heart can go at 200 bpm at most in extreme situations. I take 300 bpm
% to have a margin of 5 Hz (T=200 ms). Thus I detect peaks that are 200 ms
% apart (or 40 indexs apart) -> minpeakdistance

%% Step 3: compute thresholds for the R peaks
% First initialize the peak thresholds considering the SD of no peak areas
% and the mean of the peaks
 
thr = zeros(8528,1);
npk = zeros(8528,1);
spk = zeros(8528,1);
thr1 = zeros(8528,1);
thr2 = zeros(8528,1);

for i=1:8528
    thr(i) = mean(intA(i,:));
    nopeak = intA(i,find(intA(i,:)<thr(i)));
    npk(i) = std(nopeak);
    yespeak = intA(i,find(intA(i,:)>thr(i)));
    spk(i) = mean(yespeak);
    thr1(i)= npk(i) + 0.25*(spk(i) - npk(i));
    for q=1:10
        nopeak = intA(i,find(intA(i,:)<thr1(i)));
        npk(i) = std(nopeak);
        yespeak = intA(i,find(intA(i,:)>thr1(i)));
        spk(i) = mean(yespeak);
        thr1(i)= npk(i) + 0.25*(spk(i) - npk(i));
    end
    thr2(i) = 0.5*thr1(i);
end

%% Apply Pan-Tomkins itself to find interactive trhesholds to detect the R peaks from the processed signal

thr_i = zeros(1,8528);
pks_thr = zeros(max_pks,8528);
locs_thr = zeros(max_pks,8528);
pks_ok = zeros(8528,1);
quins_pks = zeros(8528,max_pks);
for i=1:8528
    p=0;
    npk_i = npk(i);
    spk_i = spk(i);
    thr_i(i) = npk_i + 0.25*(spk_i - npk_i);
    meandist=1000;
    for k=1:num_pks(i)
        if pks(k,i)>=thr_i(i)
            p = p+1;
            spk_i = 0.125*pks(k,i) + 0.875*spk_i;
            pks_thr(p,i) = pks(k,i);
            locs_thr(p,i) = locs(k,i);
        elseif ((pks(k,i)>=(thr_i(i)/2)) & (sum(locs_thr(:,i) ~= 0) ~= 0))
            if abs(locs(k,i) - locs_thr(p,i)) > 1.66*meandist
                p = p+1;
                spk_i = 0.125*pks(k,i) + 0.875*spk_i;
                pks_thr(p,i) = pks(k,i);
                locs_thr(p,i) = locs(k,i);
                quins_pks(i,p) = locs(k,i);
            end
        else
            npk_i = 0.125*pks(k,i) + 0.875*npk_i;
        end
        thr_i(i) = npk_i + 0.25*(spk_i - npk_i);    
        aux = locs_thr(:,i);
        meandist = (locs_thr(max(find(aux)),i) - locs_thr(1,i))/(max(find(aux))-1);
    end
    if (meandist > 75) & (meandist < 225)
        pks_ok(i) = 1;
    end
end

%% Evaluation of the detected peaks
figure(15)
aux = intA(56,:);
aux2 = locs_thr(:,56);
plot(intA(56,1:max(find(aux))))
hold on
plot(locs_thr(1:max(find(aux2)),56),pks_thr(1:max(find(aux2)),56),'o')
xlim([0 6017])
title("Incorrect peak detection of the integrated signal")
xlabel("Samples (n)")
ylabel("Amplitude")
figure(16)
aux = intA(1,:);
aux2 = locs_thr(:,1);
plot(intA(1,1:max(find(aux))))
hold on
plot(locs_thr(1:max(find(aux2)),1),pks_thr(1:max(find(aux2)),1),'o')
xlim([0 6017])
title("Correct peak detection of the integrated signal")
xlabel("Samples (n)")
ylabel("Amplitude")

%% nn or rr signals (distance between peaks)
rr = zeros(8528,(max_pks-1));
sdiff = zeros(8528,(max_pks-2)); %per calcular sdsd despres

for i=1:8528
    aux = locs_thr(:,i);
    for k=1:(max(find(aux))-1)
        rr(i,k) = locs_thr(k+1,i) - locs_thr(k,i);
    end
    aux = rr(i,:);
    for k=1:(max(find(aux))-1)
        sdiff(i,k) = rr(i,k+1) - rr(i,k);
    end
end

%% Autocorrelation of sample signals

autoAA = xcorr(AA);
autoAN = xcorr(AN);
autoAO = xcorr(AO);
autoAW = xcorr(AW);
figure(17)
subplot(4,1,1)
plot(autoAA)
subplot(4,1,2)
plot(autoAN)
subplot(4,1,3)
plot(autoAO)
subplot(4,1,4)
plot(autoAW)

%% Autocorrelation of signal A
autoA1 = zeros(1,8528);
autoA2 = zeros(1,8528);
autoA3 = zeros(1,8528);

for i=1:8528
    aux=A_re(i,:);
    autoA1(i) = mean(xcorr(A_re(i,1:max(find(aux)))));
    autoA2(i) = std(xcorr(A_re(i,1:max(find(aux)))));
    autoA3(i) = max((xcorr(A_re(i,1:max(find(aux))))));
end
autoA1(isnan(autoA1)) = 0;
autoA1(isinf(autoA1)) = 0;
autoA2(isnan(autoA2)) = 0;
autoA2(isinf(autoA2)) = 0;
autoA3(isnan(autoA3)) = 0;
autoA3(isinf(autoA3)) = 0;
autoA1 = (autoA1 - mean(autoA1))/std(autoA1);
autoA2 = (autoA2 - mean(autoA2))/std(autoA2);
autoA3 = (autoA3 - mean(autoA3))/std(autoA3);

%% Features original + autocorrelation
sdA = zeros(1,8528);
mA = zeros(1,8528);
autocorrA = zeros(1,8528);

for i=1:8528
    aux = A_re(i,:);
    sdA(i) = std(A_re(i,1:max(find(aux))));
    mA(i) = mean(A_re(i,1:max(find(aux))));
end

sdA(isnan(sdA)) = 0;
mA(isnan(mA)) = 0;
sdA(isinf(sdA)) = 0;
mA(isinf(mA)) = 0;
sdA = (sdA - mean(sdA))/std(sdA);
mA = (mA - mean(mA))/std(mA);

%% Features RR signal
sdRR = zeros(1,8528);
mRR = zeros(1,8528);
sdsd = zeros(1,8528);

for i=1:8528
    aux = rr(i,:);
    sdRR(i) = std(rr(i,1:max(find(aux))));
    mRR(i) = mean(rr(i,1:max(find(aux))));
    aux2 = sdiff(i,:);
    sdsd(i) = std(sdiff(i,1:max(find(aux2))));
end
sdRR(isnan(sdRR)) = 0;
mRR(isnan(mRR)) = 0;
sdRR(isinf(sdRR)) = 0;
mRR(isinf(mRR)) = 0;
sdsd(isinf(sdsd)) = 0;
sdsd(isnan(sdsd)) = 0;
sdRR = (sdRR - mean(sdRR))/std(sdRR);
mRR = (mRR - mean(mRR))/std(mRR);
sdsd = (sdsd - mean(sdsd))/std(sdsd);

%% Test the usefullness of features with ANOVA
reference = readtable("REFERENCE.csv");
p_mA = anova1(mA,reference.(2),'off');
p_sdA = anova1(sdA,reference.(2),'off');
p_mRR = anova1(mRR,reference.(2),'off');
p_sdRR = anova1(sdRR,reference.(2),'off');
p_pks_ok = anova1(pks_ok,reference.(2),'off');
p_sdsd = anova1(sdsd,reference.(2),'off');
p_autoA1 = anova1(autoA1,reference.(2),'off');
p_autoA2 = anova1(autoA2,reference.(2),'off'); 
p_autoA3 = anova1(autoA3,reference.(2),'off');

%% pNN50
pNN50 = zeros(8528,1);

for i=1:8528
    aux = rr(i,:);
    NN50 = 0;
    longaux=0;
    for k=1:(max(find(aux))-1)
        longaux = longaux+1;
        if (abs(rr(i,k)-rr(i,k+1))>(0.05*200))
            NN50 = NN50 + 1;
        end
    end
    pNN50(i) = 100*NN50/(longaux-1);
end

pNN50(isnan(pNN50)) = 0;
pNN50(isinf(pNN50)) = 0;
pNN50 = (pNN50-mean(pNN50))/std(pNN50);

%% pNN20
pNN20 = zeros(8528,1);

for i=1:8528
    aux = rr(i,:);
    NN20 = 0;
    longaux=0;
    for k=1:(max(find(aux))-1)
        longaux = longaux+1;
        if (abs(rr(i,k)-rr(i,k+1))>(0.02*200))
            NN20 = NN20 + 1;
        end
    end
    pNN20(i) = 100*NN20/(longaux-1);
end

pNN20(isnan(pNN20)) = 0;
pNN20(isinf(pNN20)) = 0;
pNN20 = (pNN20-mean(pNN20))/std(pNN20);

%% Evaluation of pNN20 and pNN50

p_pNN50 = anova1(pNN50,reference.(2),'off');
p_pNN20 = anova1(pNN20,reference.(2),'off');

%% pNNw (el pNN with a variable w)
pNNw = zeros(8528,20);
p_pNNw = zeros(1,20);
NNw = 0;
NN = 0;

for w=10:10:200
    for i=1:8528
        aux = rr(i,:);
        NNw = 0;
        NN = 0;
        for k=1:(max(find(aux))-1)
            NN = NN + 1;
            if (abs(rr(i,k)-rr(i,k+1))>(w*200/1000))
                NNw = NNw + 1;
            end
        end
        pNNw(i,w/10) = 100*NNw/NN;
    end
end

pNNw(isinf(pNNw)) = 0;
pNNw(isnan(pNNw)) = 0;

for i=1:20
    pNNw(:,i) = (pNNw(:,i) - mean(pNNw(:,i)))/std(pNNw(:,i));
    p_pNNw(w/10) = anova1(pNNw(:,w/10),reference.(2),'off');
end

%% PCA to check the principal components of the pNNw extracted features

[coeff,~,~,~,explained] = pca(pNNw);
PC1_pNNw = (sum(coeff(:,1).*pNNw'))';
PC2_pNNw = (sum(coeff(:,2).*pNNw'))';
PC3_pNNw = (sum(coeff(:,3).*pNNw'))';
new_pNNw = cat(2,PC1_pNNw,PC2_pNNw,PC3_pNNw);
new2_pNNw = cat(2,pNNw,new_pNNw);
[coeff2,~,~,~,explained2] = pca(new_pNNw);
[coeff3,~,~,~,explained3] = pca(new2_pNNw);

%% PCs
[pcs,scrs,~,~,pexp]=pca(pNNw);
figure(18)
pareto(pexp, ["PC1" "PC2" "PC3"],0.98)
title("Principal components to compose 99% of the variance")

%% Comprovation of PCs with ANOVA

p_PC1 = anova1(PC1_pNNw,reference.(2),'off');
p_PC2 = anova1(PC2_pNNw,reference.(2),'off');
p_PC3 = anova1(PC3_pNNw,reference.(2),'off');

%% Fraction of outliers of RR signal (+-3 times sd from the mean)
outliersRR = zeros(1,8528);

for i=1:8528
    aux = rr(i,:);
    longaux=0;
    for k=1:max(find(aux))
        longaux = longaux+1;
    end
    outliersRR(i) = 100*sum(isoutlier(rr(i,:),"mean"))/longaux;
end

outliersRR(isnan(outliersRR)) = 0;
outliersRR(isinf(outliersRR)) = 0;
outliersRR = (outliersRR - mean(outliersRR))/std(outliersRR);

%% Fraction of outliers of A signal (+-3 times sd from the mean)
outliersA = zeros(1,8528);

for i=1:8528
    aux = A_re(i,:);
    longaux=0;
    for k=1:max(find(aux))
        longaux = longaux+1;
    end
    outliersA(i) = 100*sum(isoutlier(A_re(i,:),"mean"))/longaux;
end

outliersA(isnan(outliersA)) = 0;
outliersA(isinf(outliersA)) = 0;
outliersA = (outliersA - mean(outliersA))/std(outliersA);

%% Evaluate outliers fractions with ANOVA

p_outA = anova1(outliersA,reference.(2),'off');
p_outRR = anova1(outliersRR,reference.(2),'off');

%% Poincaré map (filt sd1 i sd2) + correlation between consecutive beats length
rrn = zeros(8528,266);
rrn1 = zeros(8528,266);
tilt = zeros(1,8528);
rrn_ok = rrn;
rrn1_ok = rrn1;
sd1 = double(zeros(8528,1));
sd2 = double(zeros(8528,1));
is_ellipse = zeros(1,8528);
corrNN1 = zeros(8528,1);
corrNN2 = zeros(8528,1);
corrNN3 = zeros(8528,1);
corrlenNN = zeros(8528,1);
meanQRS = zeros(8528,10);

for i=1:8528
    aux = rr(i,:);
    rrn(i,1:(max(find(aux))-1)) = rr(i,1:(max(find(aux))-1));
    rrn1(i,1:(max(find(aux))-1)) = rr(i,2:max(find(aux)));
    
    %correlation between consecutive beats length
    corrlenNN(i) = mean(xcorr(rrn(i,1:(max(find(aux))-1)),rrn1(i,1:(max(find(aux))-1))));
    
    %poincaré
    outliers1 = isoutlier(rrn(i,1:(max(find(aux))-1)));
    outliers2 = isoutlier(rrn1(i,1:(max(find(aux))-1)));
    p=0;
    for k=1:(max(find(aux))-1)
        if (outliers1(k) == 0) & (outliers2(k) == 0)
            p = p+1;
            rrn_ok(i,p) = rrn(i,k);
            rrn1_ok(i,p) = rrn1(i,k);
        end
    end
    if p >= 5
        ellipse = fit_ellipse(rrn_ok(i,1:p),rrn1_ok(i,1:p));
        if isempty(ellipse)
            sd1(i) = 5;
            sd2(i) = 10;
        elseif isempty(ellipse.status)
            sd1(i) = ellipse.short_axis;
            sd2(i) = ellipse.long_axis;
            is_ellipse(i) = 1;
            tilt(i) = ellipse.phi;
        else
            sd1(i) = 0;
            sd2(i) = 0;
        end
    else
        sd1(i) = 5;
        sd2(i) = 10;
    end
end

% Outliers of the ellipse
for i=1:8528
    if (sd1(i) == 5) | (sd2(i) == 10)
        sd1(i) = -max(sd1);
        sd2(i) = -max(sd2);
    end
end

sd1(isnan(sd1)) = 0;
sd1(isinf(sd1)) = 0;
sd2(isnan(sd2)) = 0;
sd2(isinf(sd2)) = 0;
tilt(isnan(tilt)) = 0;
tilt(isinf(tilt)) = 0;
corrlenNN(isnan(corrlenNN))=0;
corrlenNN(isinf(corrlenNN))=0;

sd1 = (sd1 - mean(sd1))/std(sd1);
sd2 = (sd2 - mean(sd2))/std(sd2);
tilt = (tilt - mean(tilt))/std(tilt);
corrlenNN=(corrlenNN-mean(corrlenNN))/std(corrlenNN);

%%  Poincaré plot

figure(19)
aux = rr(2,:);
scatter(rrn(2,1:(max(find(aux))-1)),rrn1(2,1:(max(find(aux))-1)))
title("Poincaré plot with outliers")
xlabel("RRn(n)")
ylabel("RRn+1(n)")
figure(20)
aux = rr(1,:);
scatter(rrn(1,1:(max(find(aux))-1)),rrn1(1,1:(max(find(aux))-1)))
title("Poincaré plot")
xlabel("RRn(n)")
ylabel("RRn+1(n)")

%% Check poincaré featurws with ANOVA

p_sd1 = anova1(sd1,reference.(2),'off');
p_sd2 = anova1(sd2,reference.(2),'off');
p_tilt = anova1(tilt,reference.(2),'off');
p_corrlenNN = anova1(corrlenNN,reference.(2),'off');
p_num_pks = anova1(num_pks,reference.(2),'off');

%% bandpowers
bp1 = zeros(1,8528);
fr1 = [0.5 5];
bp2 = zeros(1,8528);
fr2 = [5 15];
bp3 = zeros(1,8528);
fr3 = [15 30];
bp4 = zeros(1,8528);
fr4 = [30 45];

for i=1:8528
    aux = A_re(i,:);
    bp1(i) = bandpower(A_re(i,1:max(find(aux))),200,fr1);
    bp2(i) = bandpower(A_re(i,1:max(find(aux))),200,fr2);
    bp3(i) = bandpower(A_re(i,1:max(find(aux))),200,fr3);
    bp4(i) = bandpower(A_re(i,1:max(find(aux))),200,fr4);    
end

%% Check bandpowers with ANOVA
p_bp1 = anova1(bp1, reference.(2),'off');
p_bp2 = anova1(bp2, reference.(2),'off');
p_bp3 = anova1(bp3, reference.(2),'off');
p_bp4 = anova1(bp4, reference.(2),'off');

%% Feature selection

% Based on the ANOVA p-values the following features are discarted: bp2,
% bp3, bp4, p_tilt, psdA, p_num_pks, p_mA, p_autoA2 and p_autoA3

% The features that remain are the following: autoA1, bp1, corrlenNN, mRR, outA, outRR, pc1_pNNw, pc2_pNNw, pc3_pNNw, pks_ok, sd1, sd2, sdRR, sdsd

SDRR=sdRR';
MRR=mRR';
outA=outliersA';
outRR=outliersRR';
TILT=tilt';
BP1=bp1';
SDSD=sdsd';
AUTOA1=autoA1';

%% Once the features are selected, the model can be selected and trained
% Data importing

reference = readtable("REFERENCE.csv", 'readvariablenames', false);
reference(:,1)=[];

%% Features 

features=reference;
features(:,"Var2")=[];
features(:,'sdRR')=array2table(SDRR);
features(:,'mRR')=array2table(MRR);
features(:,'outliersA')=array2table(outA);
features(:,'outliersRR')=array2table(outRR);
features(:,'sd1')=array2table(sd1);
features(:,'sd2')=array2table(sd2);
features(:,'pc1_pnnw')=array2table(PC1_pNNw);
features(:,'pc2_pnnw')=array2table(PC2_pNNw);
features(:,'pc3_pnnw')=array2table(PC3_pNNw);
features(:,'pks_ok')=array2table(pks_ok);
features(:,'autoA1')=array2table(AUTOA1);
features(:,'BP1')=array2table(BP1);
features(:,'corrlenNN')=array2table(corrlenNN);
features(:,'sdsd')=array2table(SDSD);

%% PCA

[pcs,scrs,~,~,pexp]=pca(table2array(features));

% Pareto plot
figure(21),
pareto(pexp)

% How many principal components for 95% of variance
pos = find(cumsum(pexp) >= 95, 1)

% Biplot
varnames={'sdRR','mRR','outlierA','outliersRR','sd1','sd2','pc1_pNNw','pc2_pNNw','pc3_pNNw','pks_ok','autoA1','bp1','corrlenNN','sdsd'};
figure(22),
biplot(pcs(:,1:2),"VarLabels",varnames)

% Heatmap
figure(23),
heatmap(abs(pcs(:,1:8)),"YDisplayLabels",varnames);
xlabel("Principal Component")

% Taking only corrlenNN
PCApreds = scrs(:,1:8);
data=array2table(PCApreds);

% Response variable
data(:,'Class')=reference(:,'Var2');

%% Train test split
cv = cvpartition(size(data,1),'HoldOut',0.3);
idx = cv.test;
train = data(~idx,:);
test  = data(idx,:);

%% Knn model
% hyperparameters already optimize

knnmodel = fitcknn(train,"Class",'NumNeighbors',14,'Distance','cityblock');
predictions_knn=predict(knnmodel,test);
labels= {'N'  'A'  'O'  '~'};
class=test.Class;
iscorrect_knn=cell2mat(predictions_knn)==cell2mat(class);
accuracy_knn=sum(iscorrect_knn)/numel(iscorrect_knn)
confusionmatrix_knn=confusionmat(class,predictions_knn);
figure(24),
confusionchart(confusionmatrix_knn,labels)

%% Svm model - fitcecoc for more than 2 classes - (VERY LOW ACCURACY!!)

svmmodel=fitcecoc(train,"Class");
predictions_svm=predict(svmmodel,test);

iscorrect_svm=cell2mat(predictions_svm)==cell2mat(class);
accuracy_svm=sum(iscorrect_svm)/numel(iscorrect_svm)
confusionmatrix_svm=confusionmat(class,predictions_svm);
figure(25),
confusionchart(confusionmatrix_svm)

%% Classification trees

treemodel=fitctree(train,"Class",'MinLeafSize',49); %optimized parameters
predictions_tree=predict(treemodel,test);

iscorrect_tree=cell2mat(predictions_tree)==cell2mat(class);
accuracy_tree=sum(iscorrect_tree)/numel(iscorrect_tree)
confusionmatrix_tree=confusionmat(class,predictions_tree);
figure(26),
confusionchart(confusionmatrix_tree,labels)

%% Separate A from B

table=reference;
table(:,'bandpower')=array2table(vectorA_re);
table(:,'Var2')=[];
table(8529:13527,'bandpower')=array2table(vectorB_re);
char='A';
column = repelem({char},8528,1);
char2='B';
column2 = repelem({char2},4999,1);
table(1:8528,'class')=column;
table(8529:13527,'class')=column2;

labels2 = {'A','B'};

%% Train test split
cv = cvpartition(size(table,1),'HoldOut',0.3);
idx = cv.test;
traindata = table(~idx,:);
testdata  = table(idx,:);

%% Knn model

knnmodel2 = fitcknn(traindata,"class",'NumNeighbors',17,'Distance','chebychev') %optimized parameters
predictions_knn2=predict(knnmodel2,testdata);
iscorrect_knn2=cell2mat(predictions_knn2)==cell2mat(testdata.class);
accuracy=sum(iscorrect_knn2)/numel(iscorrect_knn2)
confusionmatrix=confusionmat(testdata.class,predictions_knn2);
figure(27),
confusionchart(confusionmatrix,labels2)

%% Tree classifier

treemodel2=fitctree(traindata,"class",'MinLeafSize',78); %optimized parameters
predictions_tree2=predict(treemodel2,testdata);
iscorrect_tree2=cell2mat(predictions_tree2)==cell2mat(testdata.class);
accuracy=sum(iscorrect_tree2)/numel(iscorrect_tree2)
confusionmatrix=confusionmat(testdata.class,predictions_tree2);
figure(28),
confusionchart(confusionmatrix,labels2)







