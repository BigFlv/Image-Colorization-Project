close all
clear
clc

%% parametri
numVal = 20;
afisare = 1;

% parametri antrenare retea
doTraining = true;
[dataDir, ~, ~] = fileparts(mfilename('fullpath')); % incarca in dataDir calea curenta a script-ului
miniBatchSize = 30; % cate imagini ia simultan
maxEpochs = 1; % de cate ori trec setul de date prin retea
inputSize = [128, 128, 1];

% parametri pentru specificarea variantei de antrenare
% specifica varianta de antrenare
buildArchitecture = true; % !! se seteaza pe false in cazul in care au loc antrenari cu reziduuri !! 
trainGreyToLab = false;
trainGreyToLab_ch = false; chL = false; chA = false; chB = false;
trainGreyToLab_residuals = false;
trainGreyToLab_residuals_ch = false; chL_res = false; chA_res = false; chB_res = false;
trainGreyToHsv = false;
trainGreyToHsv_ch = false; chH = false; chS = false; chV = false;
trainGreyToHsv_residuals = false;
trainGreyToRgb = true;
trainGreyToRgb_residuals = false;
% specifica numele mat-ului care va contine reteaua antrenata
numeRetea = "GreyToReziduuri_50ep_PlecandDeLa_GreyToLab(cu_50ep)";

%% creare seturi de date
trainImIn = datastore(strcat(dataDir, "\train\trainSize\in"));
trainImOut = datastore(strcat(dataDir, "\train\trainSize\out"));

testImIn = datastore(strcat(dataDir, "\test\testSize\in"));
testImOut = datastore(strcat(dataDir, "\test\testSize\out"));

valImIn = datastore(strcat(dataDir, "\val\valSize\in"));
valImOut = datastore(strcat(dataDir, "\val\valSize\out"));

%% impartirea path-urilor pe diferite variante de antrenari
% Varianta 1 - path-uri pentru antrenarea retelei din gri spre RGB
if trainGreyToRgb == true
    combTrain = combine(trainImIn, trainImOut);
    combVal = combine(valImIn, valImOut);
end

% Varianta 2 - path-uri pentru antrenarea retelei din gri spre LAB
if trainGreyToLab == true
    trainImOutLabMod = datastore(strcat(dataDir, "\train\trainLab\outLab128"));
    valImOutLabMod = datastore(strcat(dataDir, "\val\valLab\outLab128"));

    combTrain = combine(trainImIn, trainImOutLabMod);
    combVal = combine(valImIn, valImOutLabMod);
end

% Varianta 3 - path-uri pentru antrenarea retelei din gri spre fiecare foaie din LAB
if trainGreyToLab_ch == true
    trainImLab = datastore(strcat(dataDir, "\train\trainLab\outLab128")); % contine imagini pe foaia Lab
    trainImOutL = datastore(strcat(dataDir, "\train\trainLab\outL")); % contine imagine pe foaia L
    trainImOutA = datastore(strcat(dataDir, "\train\trainLab\outA")); % contine imagine pe foaia A
    trainImOutB = datastore(strcat(dataDir, "\train\trainLab\outB")); % contine imagine pe foaia B

    valImLab = datastore(strcat(dataDir, "\val\valLab\outLab128"));
    valImOutL = datastore(strcat(dataDir, "\val\valLab\outL"));
    valImOutA = datastore(strcat(dataDir, "\val\valLab\outA"));
    valImOutB = datastore(strcat(dataDir, "\val\valLab\outB"));

    if chL == true
        combTrain = combine(trainImIn, trainImOutL);
        combVal = combine(valImIn, valImOutL);
    elseif chA == true
        combTrain = combine(trainImIn, trainImOutA);
        combVal = combine(valImIn, valImOutA);
    elseif chB == true
        combTrain = combine(trainImIn, trainImOutB);
        combVal = combine(valImIn, valImOutB);
    end
end

% Varianta 4 - path-uri pentru antrenarea retelei din gri spre HSV
if trainGreyToHsv == true
    trainImOutHsvMod = datastore(strcat(dataDir, "\train\trainHsv\outHsvMod"));
    valImOutHsvMod = datastore(strcat(dataDir, "\val\valHsv\outHsvMod"));

    combTrain = combine(trainImIn, trainImOutHsvMod);
    combVal = combine(valImIn, valImOutHsvMod);
end

% Varianta 5 - path-uri pentru antrenarea retelei din gri spre fiecare foaie din HSV
if trainGreyToHsv_ch == true
    trainImOutH = datastore(strcat(dataDir, "\train\trainHsv\outH")); % contine imagine pe foaia H
    trainImOutS = datastore(strcat(dataDir, "\train\trainHsv\outS")); % contine imagine pe foaia S
    trainImOutV = datastore(strcat(dataDir, "\train\trainHsv\outV")); % contine imagine pe foaia V

    valImOutH = datastore(strcat(dataDir, "\val\valHsv\outH"));
    valImOutS = datastore(strcat(dataDir, "\val\valHsv\outS"));
    valImOutV = datastore(strcat(dataDir, "\val\valHsv\outV"));

    if chH == true
        combTrain = combine(trainImIn, trainImOutH);
        combVal = combine(valImIn, valImOutH);
    elseif chS == true
        combTrain = combine(trainImIn, trainImOutS);
        combVal = combine(valImIn, valImOutS);
    elseif chV == true
        combTrain = combine(trainImIn, trainImOutV);
        combVal = combine(valImIn, valImOutV);
    end
end

% Varianta 6 - path-uri pentru antrenarea retelei din gri spre reziduurile obtinute din diferenta dintre imaginile originale si imaginile rgb generate anterior
if trainGreyToRgb_residuals == true
    trainImOut_rezRgb = datastore(strcat(dataDir,"\train\reziduuri\res_rgb\outReziduuri"));
    valImOut_rezRgb = datastore(strcat(dataDir,"\val\reziduuri\res_rgb\outReziduuri"));

    combTrain = combine(trainImIn, trainImOut_rezRgb);
    combVal = combine(valImIn, valImOut_rezRgb);
end

% Varianta 7 path-uri pentru antrenarea retelei din gri spre reziduurile obtinute din diferenta dintre imaginile lab originale si imaginile lab generate anterior
if trainGreyToLab_residuals == true
    trainImOut_rezLab = datastore(strcat(dataDir, "\train\reziduuri\res_lab\outReziduuri_Lab_pe_foi"));
    valImOut_rezLab = datastore(strcat(dataDir, "\val\reziduuri\res_lab\outReziduuri_Lab_pe_foi"));

    combTrain = combine(trainImIn, trainImOut_rezLab);
    combVal = combine(valImIn, valImOut_rezLab);
end

% Varianta 8 - path-uri pentru antrenarea retelei cu rezidurile dintre (un canal original din lab) si (un canal generat din grey to lab)
if trainGreyToLab_residuals_ch == true
    trainImOut_rez_lab_foaia_L = datastore(strcat(dataDir,"\train\reziduuri\res_lab\outReziduuri_Lab_foaia_L"));
    trainImOut_rez_lab_foaia_A = datastore(strcat(dataDir,"\train\reziduuri\res_lab\outReziduuri_Lab_foaia_A"));
    trainImOut_rez_lab_foaia_B = datastore(strcat(dataDir,"\train\reziduuri\res_lab\outReziduuri_Lab_foaia_B"));

    valImOut_rez_lab_oaia_L = datastore(strcat(dataDir,"\val\reziduuri\res_lab\outReziduuri_Lab_foaia_L"));
    valImOut_rez_lab_foaia_A = datastore(strcat(dataDir,"\val\reziduuri\res_lab\outReziduuri_Lab_foaia_A"));
    valImOut_rez_lab_foaia_B = datastore(strcat(dataDir,"\val\reziduuri\res_lab\outReziduuri_Lab_foaia_B"));
    if chL_res == true
        combTrain = combine(trainImIn, trainImOut_rez_lab_foaia_L);
        combVal = combine(valImIn, valImOut_rez_lab_oaia_L);
    elseif chA_res == true
        combTrain = combine(trainImIn, trainImOut_rez_lab_foaia_A);
        combVal = combine(valImIn, valImOut_rez_lab_foaia_A);
    elseif chB_res == true
        combTrain = combine(trainImIn, trainImOut_rez_lab_foaia_B);
        combVal = combine(valImIn, valImOut_rez_lab_foaia_B);
    end
end

% Varianta 9 - path-uri pentru antrenarea retelei cu rezidurile dintre (img hsv originale) si (img hsv generate din grey to hsv)
if trainGreyToHsv_residuals == true
    trainImOut_rezHSV = datastore(strcat(dataDir,"\train\reziduuri\res_hsv\outReziduuri_HSV"));
    valImOut_rezHSV = datastore(strcat(dataDir,"\val\reziduuri\res_hsv\outReziduuri_HSV"));
    combTrain = combine(trainImIn, trainImOut_rezHSV);
    combVal = combine(valImIn, valImOut_rezHSV);
end
dsValFull = shuffle(combVal);
dsVal = subset(dsValFull,1:numVal);
dsTrainFull = shuffle(combTrain);

% Preprocess Training and Validation Data
patchesPerImage = 12;
dsTrain = transform(dsTrainFull, ...
    @(data) extractRandomPatchV2(data,inputSize,patchesPerImage));
% Augment Training Data
dsTrain = transform(dsTrain,@(data) augmentPatchesForLowLightRecovery(data));

% afiseaza imaginea in tonuri de gri si un patch ales aleator din aceasta
if afisare == 1
    previewFull = preview(dsTrainFull);
    previewPatch = preview(dsTrain);
    figure;montage({previewFull{1,2},previewPatch{1,2}},BackgroundColor="w");
end

% dsVal contine un patch din centrul imagini
dsVal = transform(dsVal,@(data) extractRandomPatchV2(data,inputSize,patchesPerImage));

if afisare == 1
    imagePairs = read(dsTrain); % imagePairs contine patch-urile
    %x = randi(12);
    rawImage = imagePairs{1, 1}; % patch gri
    rgbPatch = imagePairs{1, 2}; % patch color
    figure;montage({rawImage,rgbPatch});
end

%% constructie arhitectura retea
if buildArchitecture == true
    load('trainedLowLightCameraPipelineNet.mat');
    netOld = netTrained;
    net = layerGraph(netOld);
    
    % strat de intrare este acelasi pt toate antrenarile
    newinputLayer = imageInputLayer(inputSize,Normalization="zerocenter");
    newinputLayer.Mean = calcMean(trainImIn);
    net = replaceLayer(net,"imageinput", newinputLayer);
    
    % modificare primul strat de convolutie
    newConvLayer = convolution2dLayer(3, 32, Padding="same", ...               
        PaddingValue="replicate",WeightsInitializer="he");
    oldLayer = netOld.Layers(2);                                             
    w = oldLayer.Weights;
    b = oldLayer.Bias;
    newConvLayer.Bias = b;
    newConvLayer.Weights = w(:, :, 1, :);
    net = replaceLayer(net,"encoderBlock1Layer1", newConvLayer);
    
    % atunci cand antrenarea se face pe o singura foaia este necesara modificarea ultimului
    % strat de comvolutie astfel incat la iesirea retelei sa existe o singura foaie
    if trainGreyToLab_ch == true || trainGreyToHsv_ch == true
        % modificare ultimul strat de convolutie
        newConvL2 = convolution2dLayer([1 1], 4, 'Name', 'StratConv' ,'Stride', [1 1],Padding="same", ...
            PaddingValue="replicate",WeightsInitializer="he");
        oldL2 = netOld.Layers(68);
        b = oldL2.Bias;
        w = oldL2.Weights;
        newConvL2.Bias = b(:, :, 1:4);
        newConvL2.Weights = w(:, :, :, 1:4);
        net = replaceLayer(net,'FinalNetworkLayer1', newConvL2);
    end
    
    % adaugare strat de redimensionare
    newlayer = resize2dLayer("scale",0.5,"Name",'mxpF');
    net = addLayers(net,newlayer);
    net = disconnectLayers(net,"FinalNetworkLayer2","ssimL1Loss");
    net = connectLayers(net,"FinalNetworkLayer2",'mxpF');
    net = connectLayers(net,"mxpF","ssimL1Loss");
end

% pentru calcularea reziduurilor nu mai este necesara modificarea arhitecturii
if trainGreyToRgb_residuals == true
    % se incarca reteaua preantrenata in etapa 1
    % load-ul obtinut din antrenarea din gri spre RGB cu 150 de epoci
    load('greyToRgb_100+50_epoci_CuValFreq250--2023-04-17-21-49-03.mat'); 
    net = layerGraph(netTrained);
elseif trainGreyToLab_residuals == true
    % se incarca reteaua preantrenata din grey spre Lab
    load('antrenareReteaGreyToLAB_100_50Epoci_FolosindLABOut128_CuValFreq250-2023-04-14-17-43-00.mat');
    net = layerGraph(netTrained);
elseif trainGreyToLab_residuals_ch == true
    % se incarca reteaua preantrenata pe un canal din grey spre Lab
    if chL_res == true
        load('antrenareFoaie_L_100_50Epoci-CuValFreq250-2023-04-15-11-32-20.mat');
    elseif chA_res == true
        load('antrenareFoaie_A_150EpociCuValFreq250-2023-03-25-09-43-18.mat'); % 150 epoci, foaia A, fara sa scad cu 128
    elseif chB_res == true
        load('antrenareFoaie_B_150+150EpociCuValFreq250-2023-03-25-23-48-15.mat'); % 150+150 epoci, foaia B
    end
    net = layerGraph(netTrained);
elseif trainGreyToHsv_residuals == true
    % se incarca reteaua preantrenata din grey spre hsv
    load('antrenareHSV_150Epoci_FolosindHSVModificat-2023-03-30-22-59-44.mat');
    net = layerGraph(netTrained);
end

%% antrenarea retelei
options = trainingOptions("adam", ...
    Plots="training-progress", ...
    MiniBatchSize=miniBatchSize, ...
    InitialLearnRate=1e-3, ...
    MaxEpochs=maxEpochs, ... 
    ValidationData=combVal,...
    ValidationFrequency=250);
if doTraining
    checkpointsDir = fullfile(dataDir,"checkpoints");
    if ~exist(checkpointsDir,"dir")
        mkdir(checkpointsDir);
    end
    options.CheckpointPath = checkpointsDir;
    netTrained = trainNetwork(dsTrain, net, options);
    modelDateTime = string(datetime("now", Format="yyyy-MM-dd-HH-mm-ss"));
    save(fullfile(dataDir, numeRetea + "-" + modelDateTime + ".mat"), "netTrained");
end

system('sleep.bat');

%% Support functions
function dataOut = extractRandomPatchV2(data,targetRAWSize,patchesPerImage)
    dataOut = cell(patchesPerImage,2);
    raw = data{1};
    rgb = data{2};
    for idx = 1:patchesPerImage
        windowGrey = randomCropWindow2d(size(raw),targetRAWSize);
        dataOut{idx,1} =  imcrop(raw,windowGrey);
        dataOut{idx,2} =  imcrop(rgb,windowGrey);
    end
end
function meanIm = calcMean(dsTrain)
    nrIm = numel(dsTrain.Files);
    meanIm = double(readimage(dsTrain,1))/nrIm;
    for i = 2 : nrIm
        meanIm = meanIm + double(readimage(dsTrain,i))/nrIm;
    end
end