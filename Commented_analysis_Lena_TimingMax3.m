%% Organise files per fish
% Get all files from current folder
matFiles = dir('Slice*analysis_matlab.mat');
matFilesNames = {matFiles.name};

% Design a regular expression that matches your naming scheme (https://regexr.com/)
fin = cellfun(@(x)regexp(x, 'fish(\d+)','tokens'), matFilesNames, 'UniformOutput', false);
names = [];
for i = 1:length(fin)
    fish = str2num(strcat(fin{i}{1}{1})); % Concatenate all the matched groups from the regex
    names(i) = fish;
end
fishList = unique(names);

clearvars i fish names fin;

% Order slices by fish
matFilesFish = [];
for individualFish = fishList
    individualFish = pad(num2str(individualFish), 2, 'left', '0');
    indexC = strfind({matFiles.name},strcat('fish', individualFish)); %Make the string match the pattern of the naming scheme you used
    matFilesFish = [matFilesFish find(not(cellfun('isempty', indexC)))];
end
matFilesOrdered = matFiles(matFilesFish); %This orders the Matlab files per fish, it helps with the indexing of the ROIs after ANTs
matFilesOrdered = rmfield(matFilesOrdered, 'isdir');
matFilesOrdered = rmfield(matFilesOrdered, 'folder');
matFilesOrdered = rmfield(matFilesOrdered, 'bytes');
matFilesOrdered = rmfield(matFilesOrdered, 'date');
matFilesOrdered = rmfield(matFilesOrdered, 'datenum');

clearvars matFilesFish indexC individualFish;

%% Decide which stimulus trains to process in directory
% Keeps ramps and not timing
rampIndex = strfind({matFilesOrdered.name}, 'Ramp');
rampInd = ~cellfun('isempty', rampIndex);
matFilesOrdered(~rampInd) = [];

clearvars rampIndex rampInd;

% Keeps A or B stimulus trains
bFind = strfind({matFilesOrdered.name}, 'B');
bIndex = ~cellfun('isempty', bFind);
matFilesOrdered(~bIndex) = [];

clearvars bFind bIndex;

%% Extract and order output from CaImAn into a matlab matrix: matFilesOrdered
% Load up the data you want from the CaImAn output (DenoisedTraces, Noise, Baseline, ROIs, Spikes and idx_components(the "good" components)) 
name = matFilesOrdered(1).name;
denoisedCalcium = load(name, 'DenoisedTraces');
denoisedCalcium = denoisedCalcium.DenoisedTraces;
denoisedCalcium(:, 1345 : end) = []; % Delete this when cropping issues no longer a problem
noise = load(name, 'Noise'); %The noise may contain the weak responses and/or the inhibition since CaImAn convolves with an exponential decay to clean up
noise = noise.Noise;
noise(:, 1345 : end) = []; % Delete this when cropping issues no longer a problem
goodComponents = load(name, 'idx_components');
goodComponents = goodComponents.idx_components + 1;
goodCalcium = denoisedCalcium(goodComponents, :);
goodNoise = noise(goodComponents, :);
matFilesOrdered(1).goodNumber = length(goodComponents);

% Loads up the ROIs and compute their centroid
Rs = load(name, 'ROIs');
Rs = Rs.ROIs; 
Rs = Rs(:, goodComponents);
corName = strrep(name, 'analysis_matlab', 'correlation');
corIm = load(corName);
corIm = corIm.Correlation_image;
dims = size(corIm);
ROI = reshape(full(Rs), dims(1), dims(2), size(Rs,2));
centroids = zeros(size(Rs, 2), 2);
for roiNb = 1:size(ROI, 3)
    progressbar([], roiNb / size(ROI, 3));
    temp = regionprops(uint16(squeeze(ROI(:,:,roiNb))) == max(max(uint16(squeeze(ROI(:, :, roiNb))))),'Centroid');    
    temp = temp.Centroid;
    centroids(roiNb, 1 : 2) = temp;
end
matFilesOrdered(1).ROIs = centroids;

clearvars noise calcium;

% Now that the initial variables are created you iterate over all the files
for i = 2 : length(matFilesOrdered)
    progressbar(i / length(matFilesOrdered), []);
    name = matFilesOrdered(i).name;
    C = load(name, 'DenoisedTraces');
    C = C.DenoisedTraces;
    if size(C, 2) > 1344
        C(:, 1345 : end) = []; % Delete this when cropping issues no longer a problem
    end
    N = load(name, 'Noise');
    N = N.Noise;
    if size(N, 2) > 1344
        N(:, 1345 : end) = []; % Delete this when cropping issues no longer a problem
    end
    Rs = load(name,'ROIs');
    Rs = Rs.ROIs;
    F = load(name,'idx_components');
    F = F.idx_components + 1;
    Rs = Rs(:, F);
    corName = strrep(name, 'analysis_matlab','correlation');
    corIm = load(corName);
    corIm = corIm.Correlation_image;
    dims = size(corIm);
    ROI = reshape(full(Rs), dims(1), dims(2), size(Rs,2));
   
    centroids = zeros(size(Rs, 2), 2);
    for roiNb = 1 : size(ROI, 3)
        progressbar([], roiNb / size(ROI, 3));
        temp = regionprops(uint16(squeeze(ROI(:, :, roiNb))) == max(max(uint16(squeeze(ROI(:, :, roiNb))))), 'Centroid');
        temp = temp.Centroid;
        centroids(roiNb, 1 : 2) = temp;
    end
    GC = C(F, :);
    GN = N(F, :);
    goodComponents = horzcat(goodComponents, F);
    goodCalcium = vertcat(goodCalcium, GC);
    goodNoise = vertcat(goodNoise, GN);
    matFilesOrdered(i).goodNumber = matFilesOrdered(i - 1).goodNumber + length(F);
    matFilesOrdered(i).ROIs = centroids;
end

clearvars temp GC C S F N name i GS GN Rs goodComponents corIm centroids corName matFilesNames matFiles ROI roiNb denoisedCalcium dims;
 
% We Z-score the data, you can choose to include the noise or not
% Your data should be in a matrix of NbNeurons x TimePoints
ZS = zscore(goodCalcium + goodNoise, 1, 2);

clearvars goodCalcium goodNoise;

%% FOR TIMING B: Flip around B scripts to remove midpoint spike/artefact
tempZS(:, 1 : 672) = ZS(:, 673 : end);
tempZS(:, 673 : 1344) = ZS(:, 1 : 672);
ZS = tempZS;

% To confirm flipped timeseries correctly
figure;
plot(mean(ZS, 1));

clearvars tempZS;

%% OPTIONAL: Detrend traces to flatten baselines to 0
figure;
plot(mean(ZS, 1));
yline(0);
xlabel('Number of frames');
ylabel('Z-scored change in fluorescence');
title("TimingA, Before detrending");
%saveas(gcf, "Fig_TimingA, Before detrending", 'svg');

detrendedZS = detrend(ZS','linear');
ZS = detrendedZS';
figure;
plot(mean(ZS,1));
yline(0);
xlabel('Number of frames');
ylabel('Z-scored change in fluorescence');
title("TimingA, After detrending");
%saveas(gcf, "Fig_TimingA, After detrending", 'svg');
%saveas(gcf, "Fig_TimingA, After detrending", 'fig');

clearvars detrendedZS;

%% Save z-scored data
save('Zscored_CaImAn_TimingA', '-v7.3');

%% Load z-scored data
load('Zscored_CaImAn_TimingA');

%% Build an index for each ROI of its corresponding fish and plane number
numbers = [0 [matFilesOrdered.goodNumber]];
idxPlane = nan(size(ZS, 2), 1);
idxFish = nan(size(ZS, 2), 1);
name = strcat(matFilesOrdered(1).name);
for i = 1 : length(matFilesOrdered)	
    name = strcat(matFilesOrdered(i).name);
    [plane, ~] = regexp(name, 'Slice(\d+)', 'tokens', 'match'); 
    plane = str2num(plane{1} {1});
    idxPlane(numbers(i) + 1: numbers(i + 1)) = plane;    
    [fish, ~] = regexp(name,'fish(\d+)','tokens'); 
    fish = str2num(fish{1} {1});
    idxFish(numbers(i) + 1: numbers(i + 1)) = fish;    
end

clearvars i fish plane name counter numbers;

%% Specify the shape of the regressors for evoked activity
loomSpike = [0, 0.419681669047195, 0.847511199979885, 1, 0.719485808200540, 0.446770218024352, 0.279790159080998, 0.178082027289836, 0.108419307355297, 0.0596454103793969, 0.0295897654706976, 0.00710019280904086, -0.00676604961064822]';
frameRate = 2;

%% FOR TIMING A: Build stimulus train based on regressor shape and timing
% List of stimulus timing for Timing A MSI
% Audio=[57,351];Loom=[15,397];-400ms=[275,522];-200ms=[100,652];-100ms=[235,482];0ms=[146,566];+100ms=[316,607];+200ms=[190,440];
burstStimulusSecs = [57, 100, 146, 190, 235, 275, 316, 351, 440, 482, 522, 566, 607, 652]; % white noise burst alone
loomStimulusSecs = [15, 100, 146, 190, 235, 275, 316, 397, 440, 482, 522, 566, 607, 652]; % checkerboard loom alone

% Overlaying the timing of Timing B MSI stimuli over the traces
stimulusTrain = zeros(2, size(ZS, 2));
for i = burstStimulusSecs
    idx = i * frameRate;
    stimulusTrain(1, idx + 9: idx - 1 + length(loomSpike) + 9) = loomSpike;
end

for i = loomStimulusSecs 
    idx = i * frameRate;
    stimulusTrain(2, idx + 9 : idx - 1 + length(loomSpike) + 9) = loomSpike;
end

% Your stimulus should look like a multicolored set of GCaMP waves
figure;
plot(stimulusTrain', 'LineWidth', 2);
hold on; 
plot(nanmean(ZS, 1), 'Color', 'k');
xlabel('Number of frames');
ylabel('Z-scored change in fluorescence');
title("TimingA, Mean Z-scored activity trace of all ROIs");
legend("Burst", "Loom");
saveas(gcf, "Fig_TimingA, Mean Z-scored activity trace of all ROIs", 'svg');
saveas(gcf, "Fig_TimingA, Mean Z-scored activity trace of all ROIs", 'fig');

clearvars loomSpike rampSpike loomStimulusSecs i idx frameRate burstStimulusSecs loomStimulusSecs;

%% FOR TIMING B: Build stimulus train based on regressor shape and timing 
% List of stimulus timing for Timing B MSI
% Audio=[15,397];Loom=[57,351];-400ms=[190,607];-200ms=[316,440];-100ms=[146,566];0ms=[235,482];+100ms=[275,652];+200ms=[100,522]
burstStimulusSecs = [15, 100, 146, 190, 235, 275, 316, 397, 440, 482, 522, 566, 607, 652]; % white noise burst alone
loomStimulusSecs = [57, 100, 146, 190, 235, 275, 316, 351, 440, 482, 522, 566, 607, 652]; % checkerboard loom alone

% Overlaying the timing of Timing B MSI stimuli over the traces
stimulusTrain = zeros(2, size(ZS, 2));
for i = burstStimulusSecs % white noise burst alone
    idx = i * frameRate;
    stimulusTrain(1, idx + 9 : idx - 1 + length(loomSpike) + 9) = loomSpike;
end

for i = loomStimulusSecs % checkerboard loom alone
    idx = i * frameRate;
    stimulusTrain(2, idx + 9 : idx - 1 + length(loomSpike) + 9) = loomSpike;
end

% Your stimulus should look like a multicolored set of GCaMP waves
figure;
plot(stimulusTrain', 'LineWidth', 2);
hold on; 
plot(mean(ZS, 1), 'Color', 'k');
xlabel('Number of frames');
ylabel('Z-scored change in fluorescence');
legend("Burst", "Loom");
title("TimingA, Mean Z-scored activity trace of all ROIs");
%saveas(gcf, "Fig_TimingA, Mean Z-scored activity trace of all ROIs", 'svg');
%saveas(gcf, "Fig_TimingA, Mean Z-scored activity trace of all ROIs", 'fig');

clearvars loomSpike rampSpike loomStimulusSecs i idx frameRate burstStimulusSecs loomStimulusSecs;

%% Calculate linear regression for each regressor (do independently using for loop)
% Perform regression on individual stimuli presentations (not combo)
linearModelStimuli = {};
allRsquaresValues = [];
for i = 1 : size(stimulusTrain, 1)
    linearModelStimuliTemp = [];
    for j = 1 : size(ZS, 1)
        mdl = fitlm(stimulusTrain(i, [1:169, 684:851])', ZS(j, [1:169, 684:851]));
        linearModelStimuliTemp(j).coef = mdl.Coefficients;
        linearModelStimuliTemp(j).rsquared = mdl.Rsquared.Adjusted;
    end
    allRsquaresValues(:, i) = [linearModelStimuliTemp.rsquared];
    linearModelStimuli{1, i} = linearModelStimuliTemp;
end

% save regression coefficient for each regressor
for i = 1 : size(linearModelStimuli, 2) 
    for j = 1 : size(ZS, 1)
        allCoefValues(j, i) = linearModelStimuli{1, i}(j).coef{2, 1}; 
    end
end

% Save all coef and rsquared values in excel
xlswrite('TimingA_allCoefValues.xlsx', allCoefValues);
xlswrite('TimingA_allRsquaresValues.xlsx', allRsquaresValues);

% Graphs the frequency distribution for rsquared values
for i = 1 : size(allRsquaresValues, 2)
    figure;
    histogram(allRsquaresValues(:, i));
    xlabel('R-squared value');
    ylabel('Number of ROIs');
    title("TimingA, Distribution of rsq for all ROIs - regressor " + i);
    saveas(gcf, "Fig_TimingA, Distribution of rsq for all ROIs - regressor " + i, 'svg');
    saveas(gcf, "Fig_TimingA, Distribution of rsq for all ROIs - regressor " + i, 'fig');
end

% Graphs the frequency distribution for regression coefficients
for i = 1 : size(allCoefValues, 2)
    figure;
    histogram(allCoefValues(:, i));
    xlabel('Regression coefficients');
    ylabel('Number of ROIs');
    title("TimingA, Distribution of coefficient for all ROIs - regressor " + i);
    saveas(gcf, "Fig_TimingA, Distribution of coefficient for all ROIs - regressor " + i, 'svg');
    saveas(gcf, "Fig_TimingA, Distribution of coefficient for all ROIs - regressor " + i, 'fig');
end

clearvars j i linearModelStimuliTemp mdl;

%% Save linear regression data 
save('LinearRegression_StimAlone_TimingA', '-v7.3');

%% Load linear regression data
load('LinearRegression_StimAlone_TimingA');

%% Applies threshold for ROI above RSQ and Coef values for ANY of the regressors (uses MAX)
rsquaresThres = 0.05; %Input
coefThres = 1; %Input

idxAboveThreshold = zeros(size(ZS, 1), 1);
for i = 1 : size(ZS, 1)
    idxAboveThreshold(i, :) = (max(allRsquaresValues(i, :)) >= rsquaresThres) && (max(allCoefValues(i, :)) >= coefThres || min(allCoefValues(i, :)) < -coefThres);
end

idxAboveThreshold = logical(idxAboveThreshold);
thresholdRsquaresValues = allRsquaresValues(idxAboveThreshold(:, 1), :);
thresholdCoefValues = allCoefValues(idxAboveThreshold(:, 1), :);
thresholdZS = ZS(idxAboveThreshold(:, 1), :);

for i = 1 : size(stimulusTrain, 1)
    figure;
    scatter(allRsquaresValues(:, i), allCoefValues(:, i),'.', 'b');
    hold on;
    scatter(thresholdRsquaresValues(:, i), thresholdCoefValues(:, i),'.', 'r');
    axis([-0.1 1 -6 12]);
    xlabel('R-sqaured values');
    ylabel('Regression coefficients');
    title("TimingA, All ROIs above thresholds " + rsquaresThres + " and " + coefThres + " for regressor " + i);
    saveas(gcf, "Fig_TimingA, Rsq versus coefficient for all ROIs - regressor " + i, 'svg');
    saveas(gcf, "Fig_TimingA, Rsq versus coefficient for all ROIs - regressor " + i, 'png');
    saveas(gcf, "Fig_TimingA, Rsq versus coefficient for  all ROIs - regressor " + i, 'fig');
end

clearvars i rsquaresThres coefThres thresholdCoefValues thresholdRsquaresValues;

%% Generates index for fish and planes for new thresholded ROIs
idxFishZS = idxFish(idxAboveThreshold); 
idxPlaneZS = idxPlane(idxAboveThreshold);

clearvars allRsquaresValues allCoefValues;

%% FOR TIMING A: Inputs required to calculate max responses
% Input the times in seconds when stimuli come on
burstMaskSecsR1 = [124, 135];
loomMaskSecsR1 = [40, 51];
neg400MaskSecsR1 = [560, 571];
neg200MaskSecsR1 = [210, 221];
neg100MaskSecsR1 = [480, 491];
zeroMaskSecsR1 = [302, 313];
pos100MaskSecsR1 = [642, 653];
pos200MaskSecsR1 = [390, 401];

burstMaskSecsR2 = [712, 723];
loomMaskSecsR2 = [804, 815];
neg400MaskSecsR2 = [1054, 1065];
neg200MaskSecsR2 = [1314, 1325];
neg100MaskSecsR2 = [974, 985];
zeroMaskSecsR2 = [1142, 1153];
pos100MaskSecsR2 = [1224, 1235];
pos200MaskSecsR2 = [890, 901];

%% FOR TIMING B: Inputs required to calculate max responses
% Input the times in seconds when stimuli come on
burstMaskSecsR1 = [40, 51];
loomMaskSecsR1 = [124, 135];
neg400MaskSecsR1 = [390, 401];
neg200MaskSecsR1 = [642, 653];
neg100MaskSecsR1 = [302, 313];
zeroMaskSecsR1 = [480, 491];
pos100MaskSecsR1 = [560, 571];
pos200MaskSecsR1 = [210, 221];

burstMaskSecsR2 = [804, 815];
loomMaskSecsR2 = [712, 723];
neg400MaskSecsR2 = [1224, 1235];
neg200MaskSecsR2 = [890, 901];
neg100MaskSecsR2 = [1142, 1153];
zeroMaskSecsR2 = [974, 985];
pos100MaskSecsR2 = [1314, 1325];
pos200MaskSecsR2 = [1054, 1065];

%% Extracts the ZS from timeseries during stimulus presentations
% Specify frames from ZS when stimuli where shown (for first block, R1)
msiBurstAloneR1 = thresholdZS(:, burstMaskSecsR1(1): burstMaskSecsR1(2));
msiLoomAloneR1 = thresholdZS(:, loomMaskSecsR1(1) : loomMaskSecsR1(2));
msiNeg400R1 = thresholdZS(:, neg400MaskSecsR1(1) : neg400MaskSecsR1(2));
msiNeg200R1 = thresholdZS(:, neg200MaskSecsR1(1) : neg200MaskSecsR1(2));
msiNeg100R1 = thresholdZS(:, neg100MaskSecsR1(1) : neg100MaskSecsR1(2));
msiZeroR1 = thresholdZS(:, zeroMaskSecsR1(1) : zeroMaskSecsR1(2));
msiPos100R1 = thresholdZS(:, pos100MaskSecsR1(1) : pos100MaskSecsR1(2));
msiPos200R1 = thresholdZS(:, pos200MaskSecsR1(1) : pos200MaskSecsR1(2));

% Specify frames from ZS when stimuli where shown (for second block, R2)
msiBurstAloneR2 = thresholdZS(:, burstMaskSecsR2(1): burstMaskSecsR2(2));
msiLoomAloneR2 = thresholdZS(:, loomMaskSecsR2(1) : loomMaskSecsR2(2));
msiNeg400R2 = thresholdZS(:, neg400MaskSecsR2(1) : neg400MaskSecsR2(2));
msiNeg200R2 = thresholdZS(:, neg200MaskSecsR2(1) : neg200MaskSecsR2(2));
msiNeg100R2 = thresholdZS(:, neg100MaskSecsR2(1) : neg100MaskSecsR2(2));
msiZeroR2 = thresholdZS(:, zeroMaskSecsR2(1) : zeroMaskSecsR2(2));
msiPos100R2 = thresholdZS(:, pos100MaskSecsR2(1) : pos100MaskSecsR2(2));
msiPos200R2 = thresholdZS(:, pos200MaskSecsR2(1) : pos200MaskSecsR2(2));

%% Calculate maximum ZS value during stimulus presentation masks
multipleSD = 1; %INPUT
maxValues = {};
for i = 1 : length(thresholdZS)
    maxValues(i).maxBurstAloneR1 = max(msiBurstAloneR1(i, :));
    maxValues(i).maxBurstAloneR2 = max(msiBurstAloneR2(i, :));
    maxValues(i).maxLoomAloneR1 = max(msiLoomAloneR1(i, :));
    maxValues(i).maxLoomAloneR2 = max(msiLoomAloneR2(i, :));
    maxValues(i).avgMaxBurstR1R2 = ([maxValues(i).maxBurstAloneR1] + [maxValues(i).maxBurstAloneR2]) / 2;
    maxValues(i).avgMaxLoomR1R2 = ([maxValues(i).maxLoomAloneR1] + [maxValues(i).maxLoomAloneR2]) / 2;
    maxValues(i).stdMaxBurstR1R2 = (abs([maxValues(i).maxBurstAloneR1] - [maxValues(i).maxBurstAloneR2]) / sqrt(2)) * multipleSD;
    maxValues(i).stdMaxLoomR1R2 = (abs([maxValues(i).maxLoomAloneR1] - [maxValues(i).maxLoomAloneR2]) / sqrt(2)) * multipleSD;
    maxValues(i).maxNeg400R1 = max(msiNeg400R1(i, :));
    maxValues(i).maxNeg400R2 = max(msiNeg400R2(i, :));
    maxValues(i).maxNeg200R1 = max(msiNeg200R1(i, :));
    maxValues(i).maxNeg200R2 = max(msiNeg200R2(i, :));
    maxValues(i).maxNeg100R1 = max(msiNeg100R1(i, :));
    maxValues(i).maxNeg100R2 = max(msiNeg100R2(i, :));
    maxValues(i).maxZeroR1 = max(msiZeroR1(i, :));
    maxValues(i).maxZeroR2 = max(msiZeroR2(i, :));
    maxValues(i).maxPos100R1 = max(msiPos100R1(i, :));
    maxValues(i).maxPos100R2 = max(msiPos100R2(i, :));
    maxValues(i).maxPos200R1 = max(msiPos200R1(i, :));
    maxValues(i).maxPos200R2 = max(msiPos200R2(i, :));
end

clearvars i multipleSD;
clearvars burstMaskSecsR1 loomMaskSecsR1 neg400MaskSecsR1 neg200MaskSecsR1 neg100MaskSecsR1 zeroMaskSecsR1 pos100MaskSecsR1 pos200MaskSecsR1;
clearvars burstMaskSecsR2 loomMaskSecsR2 neg400MaskSecsR2 neg200MaskSecsR2 neg100MaskSecsR2 zeroMaskSecsR2 pos100MaskSecsR2 pos200MaskSecsR2;
clearvars msiBurstAloneR1 msiLoomAloneR1 msiNeg400R1 msiNeg200R1 msiNeg100R1 msiZeroR1 msiPos100R1 msiPos200R1;
clearvars msiBurstAloneR2 msiLoomAloneR2 msiNeg400R2 msiNeg200R2 msiNeg100R2 msiZeroR2 msiPos100R2 msiPos200R2;

%% Define which ROIs are supralinear and sublinear on a scale of 1-10
% Returns 1 if ROI is greater than sum of unisensory stimuli + 1SD
isSupraLinear = zeros(length(thresholdZS), 10);
for i = 1 : length(thresholdZS)
    if (((maxValues(i).avgMaxLoomR1R2 + maxValues(i).stdMaxLoomR1R2) < maxValues(i).maxNeg200R1) && (maxValues(i).avgMaxBurstR1R2 + maxValues(i).stdMaxBurstR1R2) < maxValues(i).maxNeg200R1)
        isSupraLinear(i, 1) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 + maxValues(i).stdMaxLoomR1R2) < maxValues(i).maxNeg200R2) && (maxValues(i).avgMaxBurstR1R2 + maxValues(i).stdMaxBurstR1R2) < maxValues(i).maxNeg200R2)
        isSupraLinear(i, 2) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 + maxValues(i).stdMaxLoomR1R2) < maxValues(i).maxNeg100R1) && (maxValues(i).avgMaxBurstR1R2 + maxValues(i).stdMaxBurstR1R2) < maxValues(i).maxNeg100R1)
        isSupraLinear(i, 3) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 + maxValues(i).stdMaxLoomR1R2) < maxValues(i).maxNeg100R2) && (maxValues(i).avgMaxBurstR1R2 + maxValues(i).stdMaxBurstR1R2) < maxValues(i).maxNeg100R2)
        isSupraLinear(i, 4) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 + maxValues(i).stdMaxLoomR1R2) < maxValues(i).maxZeroR1) && (maxValues(i).avgMaxBurstR1R2 + maxValues(i).stdMaxBurstR1R2) < maxValues(i).maxZeroR1)
        isSupraLinear(i, 5) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 + maxValues(i).stdMaxLoomR1R2) < maxValues(i).maxZeroR2) && (maxValues(i).avgMaxBurstR1R2 + maxValues(i).stdMaxBurstR1R2) < maxValues(i).maxZeroR2)
        isSupraLinear(i, 6) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 + maxValues(i).stdMaxLoomR1R2) < maxValues(i).maxPos100R1) && (maxValues(i).avgMaxBurstR1R2 + maxValues(i).stdMaxBurstR1R2) < maxValues(i).maxPos100R1)
        isSupraLinear(i, 7) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 + maxValues(i).stdMaxLoomR1R2) < maxValues(i).maxPos100R2) && (maxValues(i).avgMaxBurstR1R2 + maxValues(i).stdMaxBurstR1R2) < maxValues(i).maxPos100R2)
        isSupraLinear(i, 8) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 + maxValues(i).stdMaxLoomR1R2) < maxValues(i).maxPos200R1) && (maxValues(i).avgMaxBurstR1R2 + maxValues(i).stdMaxBurstR1R2) < maxValues(i).maxPos200R1)
        isSupraLinear(i, 9) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 + maxValues(i).stdMaxLoomR1R2) < maxValues(i).maxPos200R2) && (maxValues(i).avgMaxBurstR1R2 + maxValues(i).stdMaxBurstR1R2) < maxValues(i).maxPos200R2)
        isSupraLinear(i, 10) = 1;
    end
end

% Returns 1 if ROI is less than sum of unisensory stimuli - 1SD
isSubLinear = zeros(length(thresholdZS), 10);
for i = 1 : length(thresholdZS)
        if (((maxValues(i).avgMaxLoomR1R2 - maxValues(i).stdMaxLoomR1R2) > maxValues(i).maxNeg200R1) && (maxValues(i).avgMaxBurstR1R2 - maxValues(i).stdMaxBurstR1R2) > maxValues(i).maxNeg200R1)
        isSubLinear(i, 1) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 - maxValues(i).stdMaxLoomR1R2) > maxValues(i).maxNeg200R2) && (maxValues(i).avgMaxBurstR1R2 - maxValues(i).stdMaxBurstR1R2) > maxValues(i).maxNeg200R2)
        isSubLinear(i, 2) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 - maxValues(i).stdMaxLoomR1R2) > maxValues(i).maxNeg100R1) && (maxValues(i).avgMaxBurstR1R2 - maxValues(i).stdMaxBurstR1R2) > maxValues(i).maxNeg100R1)
        isSubLinear(i, 3) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 - maxValues(i).stdMaxLoomR1R2) > maxValues(i).maxNeg100R2) && (maxValues(i).avgMaxBurstR1R2 - maxValues(i).stdMaxBurstR1R2) > maxValues(i).maxNeg100R2)
        isSubLinear(i, 4) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 - maxValues(i).stdMaxLoomR1R2) > maxValues(i).maxZeroR1) && (maxValues(i).avgMaxBurstR1R2 - maxValues(i).stdMaxBurstR1R2) > maxValues(i).maxZeroR1)
        isSubLinear(i, 5) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 - maxValues(i).stdMaxLoomR1R2) > maxValues(i).maxZeroR2) && (maxValues(i).avgMaxBurstR1R2 - maxValues(i).stdMaxBurstR1R2) > maxValues(i).maxZeroR2)
        isSubLinear(i, 6) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 - maxValues(i).stdMaxLoomR1R2) > maxValues(i).maxPos100R1) && (maxValues(i).avgMaxBurstR1R2 - maxValues(i).stdMaxBurstR1R2) > maxValues(i).maxPos100R1)
        isSubLinear(i, 7) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 - maxValues(i).stdMaxLoomR1R2) > maxValues(i).maxPos100R2) && (maxValues(i).avgMaxBurstR1R2 - maxValues(i).stdMaxBurstR1R2) > maxValues(i).maxPos100R2)
        isSubLinear(i, 8) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 - maxValues(i).stdMaxLoomR1R2) > maxValues(i).maxPos200R1) && (maxValues(i).avgMaxBurstR1R2 - maxValues(i).stdMaxBurstR1R2) > maxValues(i).maxPos200R1)
        isSubLinear(i, 9) = 1;
    end
    if (((maxValues(i).avgMaxLoomR1R2 - maxValues(i).stdMaxLoomR1R2) > maxValues(i).maxPos200R2) && (maxValues(i).avgMaxBurstR1R2 - maxValues(i).stdMaxBurstR1R2) > maxValues(i).maxPos200R2)
        isSubLinear(i, 10) = 1;
    end
end

% For each ROI, sums the number of supra or sub-linear responses to the 10
% stimuli (-200ms to +200ms)
sumSupraLinear = sum(isSupraLinear, 2);
sumSubLinear = sum(isSubLinear, 2); 
edges = (1 : 1 : 10);
figure;
histogram(sumSupraLinear, edges, 'FaceColor', 'g', 'FaceAlpha', 0.25);
hold on;
histogram(sumSubLinear, edges, 'FaceColor', 'r', 'FaceAlpha', 0.25);
xlabel('Number of responses to the 10 stimuli');
ylabel('Number of ROIs');
legend('Supralinear', 'Sublinear');
title("TimingA, Number of times ROI responds to 10 stimuli");
saveas(gcf, "Fig_TimingA, Number of times ROI responds to 10 stimuli", 'fig');
saveas(gcf, "Fig_TimingA, Number of times ROI responds to 10 stimuli", 'png');
saveas(gcf, "Fig_TimingA, Number of times ROI responds to 10 stimuli", 'svg');

%clearvars i edges isSupraLinear isSubLinear;

%% Extracts which auditory responsive ROIs are supra/sub-linear
% Input cutoffs
cutoffSupra = 6;
cutoffSub = 6; 

idxSupraLinear = sumSupraLinear >= cutoffSupra;
idxSubLinear = sumSubLinear >= cutoffSub;
idxNotLinear = zeros(length(thresholdZS), 1);
for i = 1 : length(thresholdZS)
    idxNotLinear(i, :) =  max(idxSupraLinear(i, :),  idxSubLinear(i, :));
end

supraLinearZS = thresholdZS(idxSupraLinear(:, 1), :);
subLinearZS = thresholdZS(idxSubLinear(:, 1), :);

% Generates index for fish and planes for new thresholded ROIs
idxNotLinear = logical(idxNotLinear);
idxFishMax = idxFishZS(idxNotLinear);
idxPlaneMax = idxPlaneZS(idxNotLinear);

clearvars i sumSuparLinear sumSubLinear;

%% Finds nonlinear ROIs
% Calculates number SupraLinear ROIs per fish
idxFishSupra = idxFishZS(idxSupraLinear);
numSupraROIperFish = zeros(length(fishList), 1);
for i = 1 : length(fishList) 
    numSupraROIperFish(i, 1) = size(supraLinearZS(idxFishSupra == fishList(i)), 1);
end

% Calculates number SubLinear ROIs per fish
idxFishSub = idxFishZS(idxSubLinear);
numSubROIperFish = zeros(length(fishList), 1);
for i = 1 : length(fishList) 
    numSubROIperFish(i, 1) = size(subLinearZS(idxFishSub == fishList(i)), 1);
end

%% FOR TIMING A: Plots summary of nonlinear ROIs
figure;
subplot(3, 1, 1);
bar(numSupraROIperFish(:, 1));
title("TimingA, Number of supralinear ROIs per fish at cutoff " + cutoffSupra);
xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]);
xticklabels({'WT', 'h', 'FMR', 'h', 'WT', 'h', 'h', 'h', 'h', 'WT', 'WT', 'h', 'h', 'WT', 'WT', 'h', 'h', 'h', 'h', 'WT', 'h', 'h', 'h', 'FMR'});
subplot(3, 1, 2);
imagesc(supraLinearZS, [-3, 8]);
title("TimingA, Z-scored activity traces of " + size(supraLinearZS, 1) + " supralinear ROIs");
subplot(3, 1, 3);
plot(mean(supraLinearZS, 1));
xlim([0 1344]);
xticks([39 123 209 301 389 479 559 641 711 803 889 973 1053 1141 1223 1313]);
xticklabels({'V', 'A', '-200', '0', '+200', '-100', '-400', '+100', 'A', 'V', '+200', '-100', '-400', '0', '+100', '-200'});
title("TimingA, Mean z-scored activity trace of " + size(supraLinearZS, 1) + " supralinear ROIs");
saveas(gcf, "Fig_TimingA, SupraLinear ROIs at cutoff " + cutoffSupra, 'svg');
saveas(gcf, "Fig_TimingA, SupraLinear ROIs at cutoff " + cutoffSupra, 'png');
saveas(gcf, "Fig_TimingA, SupraLinear ROIs at cutoff " + cutoffSupra, 'fig');

figure;
title("TimingA, Sublinear ROIs at cutoff " + cutoffSub);
subplot(3, 1, 1);
bar(numSubROIperFish(:, 1));
title("TimingA, Number of sublinear ROIs per fish at cutoff " + cutoffSub);
xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]);
xticklabels({'WT', 'h', 'FMR', 'h', 'WT', 'h', 'h', 'h', 'h', 'WT', 'WT', 'h', 'h', 'WT', 'WT', 'h', 'h', 'h', 'h', 'WT', 'h', 'h', 'h', 'FMR'});
subplot(3, 1, 2);
imagesc(subLinearZS, [-3, 8]);
title("TimingA, Z-scored activity traces of " + size(subLinearZS, 1) + " sublinear ROIs");
subplot(3, 1, 3);
plot(mean(subLinearZS));
xlim([0 1344]);
xticks([39 123 209 301 389 479 559 641 711 803 889 973 1053 1141 1223 1313]);
xticklabels({'V', 'A', '-200', '0', '+200', '-100', '-400', '+100', 'A', 'V', '+200', '-100', '-400', '0', '+100', '-200'});
title("TimingA, Mean z-scored activity trace of " + size(subLinearZS, 1) + " sublinear ROIs");
saveas(gcf, "Fig_TimingA, SubLinear ROIs at cutoff " + cutoffSub, 'svg');
saveas(gcf, "Fig_TimingA, SubLinear ROIs at cutoff " + cutoffSub, 'png');
saveas(gcf, "Fig_TimingA, SubLinear ROIs at cutoff " + cutoffSub, 'fig');

%% FOR TIMING B: Plots summary of nonlinear ROIs
figure;
subplot(3, 1, 1);
bar(numSupraROIperFish(:, 1));
title("TimingA, Number of supralinear ROIs per fish at cutoff " + cutoffSupra);
xticklabels({'het', 'het', 'het', 'het', 'het', 'het', 'FMR', 'het', 'het', 'het', 'WT', 'FMR', 'FMR', 'het'});
subplot(3, 1, 2);
imagesc(supraLinearZS, [-3, 8]);
title("TimingA, Z-scored activity traces of " + size(supraLinearZS, 1) + " supralinear ROIs");
subplot(3, 1, 3);
plot(mean(supraLinearZS, 1));
xlim([0 1344]);
xticks([39 123 209 301 389 479 559 641 711 803 889 973 1053 1141 1223 1313]);
xticklabels({'A', 'V', '+200', '-100', '-400', '0', '+100', '-200', 'V', 'A', '-200', '0', '+200', '-100', '-400', '+200'});
title("TimingA, Mean z-scored activity trace of " + size(supraLinearZS, 1) + " supralinear ROIs");
saveas(gcf, "Fig_TimingA, SupraLinear ROIs at cutoff " + cutoffSupra, 'svg');
saveas(gcf, "Fig_TimingA, SupraLinear ROIs at cutoff " + cutoffSupra, 'png');
saveas(gcf, "Fig_TimingA, SupraLinear ROIs at cutoff " + cutoffSupra, 'fig');

figure;
title("TimingA, Sublinear ROIs at cutoff " + cutoffSub);
subplot(3, 1, 1);
bar(numSubROIperFish(:, 1));
title("TimingA, Number of sublinear ROIs per fish at cutoff " + cutoffSub);
xticklabels({'het', 'het', 'het', 'het', 'het', 'het', 'FMR', 'het', 'het', 'het', 'WT', 'FMR', 'FMR', 'het'});
subplot(3, 1, 2);
imagesc(subLinearZS, [-3, 8]);
title("TimingA, Z-scored activity traces of " + size(subLinearZS, 1) + " sublinear ROIs");
subplot(3, 1, 3);
plot(mean(subLinearZS));
xlim([0 1344]);
xticks([39 123 209 301 389 479 559 641 711 803 889 973 1053 1141 1223 1313]);
xticklabels({'A', 'V', '+200', '-100', '-400', '0', '+100', '-200', 'V', 'A', '-200', '0', '+200', '-100', '-400', '+200'});
title("TimingA, Mean z-scored activity trace of " + size(subLinearZS, 1) + " sublinear ROIs");
saveas(gcf, "Fig_TimingA, SubLinear ROIs at cutoff " + cutoffSub, 'svg');
saveas(gcf, "Fig_TimingA, SubLinear ROIs at cutoff " + cutoffSub, 'png');
saveas(gcf, "Fig_TimingA, SubLinear ROIs at cutoff " + cutoffSub, 'fig');

%%
%Calculates the number of ROIs detected per fish
for i = 1 : max(idxFish) %Last ID of fish!
    fishX=sum(idxFish==i);
    numTotalROIperFish(i, :)=sum(fishX);
end

% Save all counts per fish information
xlswrite('TimingA_numSupraROIperFish.xlsx', numSupraROIperFish);
xlswrite('TimingA_numSubROIperFish.xlsx', numSubROIperFish);
xlswrite('TimingA_numTotalROIperFish.xlsx', numTotalROIperFish);

clearvars i fishX;

%% Plot of where supralinear and sublinear ROIs are located in xy
% Concatinates the coordinates of ALL ROIs
centroids = [];
for i = 1 : length(matFilesOrdered)
    centroids = cat(1, centroids, matFilesOrdered(i).ROIs);
end

centroidsThreshold = centroids(idxAboveThreshold, :);
centroidsSupra = centroidsThreshold(idxSupraLinear, :);
centroidsSub = centroidsThreshold(idxSubLinear, :);
idxRandom = randsample(size(thresholdZS, 1), 20000);

figure;
scatter(centroidsThreshold(idxRandom, 2), centroidsThreshold(idxRandom, 1), 8, [0.8, 0.8, 0.8], 'filled');
hold on
scatter(centroidsSupra(:, 2), centroidsSupra(:, 1), 8, 'filled');
scatter(centroidsSub(:, 2), centroidsSub(:, 1), 8, 'filled');
xlim([0,600]);ylim([0,600]);
title("TimingA, Location of " + size(idxFishMax, 1) + " nonlinear ROIs at cutoff " + cutoffSupra + " n " + cutoffSub);
legend('Random', 'Supralinear', 'Sublinear');
saveas(gcf, "Fig_TimingA, Location of nonlinear ROIs at cutoff " + cutoffSupra + "n" + cutoffSub, 'svg');
saveas(gcf, "Fig_TimingA, Location of nonlinear ROIs at cutoff " + cutoffSupra + "n" + cutoffSub, 'png');
saveas(gcf, "Fig_TimingA, Location of nonlinear ROIs at cutoff " + cutoffSupra + "n" + cutoffSub, 'fig');

clearvars i c;