close all
clear
clc

% Script destinat manipularii seturilor de date

%% parametri 
[dataDir, ~, ~] = fileparts(mfilename('fullpath')); % incarca in dataDir calea curenta a script-ului
resizeImg = false; newSize = [128 128];
rgbToGrey = false;
rgbToLab = false;
rgbToHsv = false;
storeGenImg = false;
target_residues = false;
target_residues_ch = false;
storesLclImgGen = true;
storesLclImgGen_ch = false;

% seturile de date originale nemodificate
trainImOrg = datastore(strcat(dataDir, "\train\out"));
testImOrg = datastore(strcat(dataDir, "\test\out"));
valImOrg = datastore(strcat(dataDir, "\val\out"));

% seturi de date rescalate
trainImIn = datastore(strcat(dataDir, "\train\trainSize\in"));
trainImOut = datastore(strcat(dataDir, "\train\trainSize\out"));

testImIn = datastore(strcat(dataDir, "\test\testSize\in"));
testImOut = datastore(strcat(dataDir, "\test\testSize\out"));

valImIn = datastore(strcat(dataDir, "\val\valSize\in"));
valImOut = datastore(strcat(dataDir, "\val\valSize\out"));

%% apeluri de functii
% Redimensionare imagini originale si le pun in X_Size\out
if resizeImg == true
    resizeImds(trainImOrg, dataDir + "\train\trainSize\out", newSize);
    resizeImds(testImOrg, dataDir + "\test\testSize\out", newSize);
    resizeImds(valImOrg, dataDir + "\val\valSize\out", newSize);
end

% Transformare imagini redimensionate in gri si le pun in X\in 
if rgbToGrey == true
    fromRgbToGrey(trainImOut, dataDir + "\train\trainSize\in");
    fromRgbToGrey(testImOut, dataDir + "\test\testSize\in");
    fromRgbToGrey(valImOut, dataDir + "\val\valSize\in");
end

% Transformare imagini din rgb in lab si stocheaza local fiecare canal
if rgbToLab == true
    convertionFromRgbToLab(trainImOut, dataDir + "\train\trainLab\Out");
    convertionFromRgbToLab(testImOut, dataDir + "\test\testLab\Out");
    convertionFromRgbToLab(valImOut, dataDir + "\val\valLab\Out");
end

% Transformare imagini din rgb in hsv si stocheaza local fiecare canal
if rgbToHsv == true
    conversionFromRgbToHsv(trainImOut, dataDir + "\train\trainHsv\Out");
    conversionFromRgbToHsv(testImOut, dataDir + "\test\testHsv\Out");
    conversionFromRgbToHsv(valImOut, dataDir + "\val\valHsv\Out");
end

% Stocarea locala a imaginilor generate
if storesLclImgGen == true
    % format Lab
    load('antrenareReteaGreyToLAB_100_50Epoci_FolosindLABOut128_CuValFreq250-2023-04-14-17-43-00.mat');
    netOld_LAB_Mod = netTrained;
    % format HSV
    load('antrenareHSV_150Epoci_FolosindHSVModificat-2023-03-30-22-59-44.mat');
    netOld_HSV_Mod = netTrained;

    storesLocallyImagesGenerated(trainImIn, dataDir, "\train\stocareLabGen", netOld_LAB_Mod);
    storesLocallyImagesGenerated(testImIn, dataDir, "\test\stocareLabGen", netOld_LAB_Mod);
    storesLocallyImagesGenerated(valImIn, dataDir, "\val\stocareLabGen", netOld_LAB_Mod);

    storesLocallyImagesGenerated(trainImIn, dataDir, "\train\stocareHSVGen", netOld_HSV_Mod);
    storesLocallyImagesGenerated(testImIn, dataDir, "\test\stocareHSVGen", netOld_HSV_Mod);
    storesLocallyImagesGenerated(valImIn, dataDir, "\val\stocareHSVGen", netOld_HSV_Mod);
end

% Formarea setului de date target de tip reziduu
if target_residues == true
    % imagini lab nemodificate
    trainImOut_Nemod_Lab = datastore(strcat(dataDir, "\train\trainLab\outLab128"));
    testImOut_Nemod_Lab = datastore(strcat(dataDir, "\test\testLab\outLab128"));
    valImOut_Nemod_Lab = datastore(strcat(dataDir, "\val\valLab\outLab128"));
    % imagini hsv nemodificate
    trainImOut_Nemod_HSV = datastore(strcat(dataDir, "\train\trainHsv\outHsvMod"));
    testImOut_Nemod_HSV = datastore(strcat(dataDir, "\test\testHsv\outHsvMod"));
    valImOut_Nemod_HSV = datastore(strcat(dataDir, "\val\valHsv\outHsvMod"));
    
    % imagini lab generate
    trainImOut_Gen = datastore(strcat(dataDir, "\train\stocareLabGen"));
    testImOut_Gen = datastore(strcat(dataDir, "\test\stocareLabGen"));
    valImOut_Gen = datastore(strcat(dataDir, "\val\stocareLabGen"));
    % imagini hsv generate
    trainImOut_Gen_HSV = datastore(strcat(dataDir, "\train\stocareHSVGen"));
    testImOut_Gen_HSV = datastore(strcat(dataDir, "\test\stocareHSVGen"));
    valImOut_Gen_HSV = datastore(strcat(dataDir, "\val\stocareHSVGen"));

    targetFormation_residues(trainImOut_Nemod_Lab, trainImOut_Gen, dataDir, "\train\reziduuri\res_lab\outReziduuri_Lab", "Lab");
    targetFormation_residues(testImOut_Nemod_Lab, testImOut_Gen, dataDir, "\test\reziduuri\res_lab\outReziduuri_Lab", "Lab");
    targetFormation_residues(valImOut_Nemod_Lab, valImOut_Gen, dataDir, "\val\reziduuri\res_lab\outReziduuri_Lab", "Lab");

    targetFormation_residues(trainImOut_Nemod_HSV, trainImOut_Gen_HSV, dataDir, "\train\reziduuri\res_hsv\outReziduri_HSV", "HSV");
    targetFormation_residues(testImOut_Nemod_HSV, testImOut_Gen_HSV, dataDir, "\test\reziduuri\res_hsv\outReziduri_HSV", "HSV"); 
    targetFormation_residues(valImOut_Nemod_HSV, valImOut_Gen_HSV, dataDir, "\val\reziduuri\res_hsv\outReziduri_HSV", "HSV");

 end

function resizeImds(imdsBefore, pathAfter, newSize)
    % functie utilizata pentru redimensionarea imaginilor dintr-un datastore
    % imdsBefore - este un datastore cu imagini care necesita redimensionate, 
    % pathAfter - este o cale pt stocare locala a seturilor de imagini redimensionate, 
    % newSize - este noua dimensiune a setului de date
    for i = 1 : numel(imdsBefore)
        %read the file name
        oldFileName=imdsBefore.Files{i};
        
        %resize the image
        im = imread(oldFileName);    
        imResized= imresize(im, newSize);
    
        %scale the image
        imScale = double(imResized) / 255;
    
        % create the new name
        [~,fName, fileExt] = fileparts(oldFileName);
        newFileName = strcat(pathAfter, '\', fName, fileExt);
    
        % write the resize image
        imwrite(imScale, newFileName);
    end
end

function fromRgbToGrey(imdsRgb, pathGrey)
    % functie utilizata pentru conversia din RGB in tonuri de gri
    % primul parametru este un datastore cu imagini color, iar al doilea este o cale pt stocare locala a imaginilor gri
    for i = 1 : numel(imdsRgb.Files)
        im = imread(imdsRgb.Files{i});
        imGrey = rgb2gray(im);

        [~, name, ext] = fileparts(imdsRgb.Files{i});
        newFileName = strcat(pathGrey, '\', name, ext);
        imwrite(imGrey, newFileName);
    end
end

function convertionFromRgbToLab(imdsRgb, pathLab)
    % Functie utilizata la conversia imaginilor din rgb in lab si stocarea lor locala. Canalele color space-ul A si B au valori in intervalul
    % (-128, 127) asadar este necesara adunarea canalelor A si B cu 128 inainte de a le stoca local pt a evita pierderea pixelilor negativi.
    for i = 1 : numel(imdsRgb.Files)
        oldFileName = imdsRgb.Files{i};
           
        im = imread(oldFileName);
        im = rgb2lab(im);
        
        %am adunat cu 128 la A si la B deoarece stocarea locala transforma imaginile in uint8 si asa se pierd pixelii negativi
        imL = uint8(im(:, :, 1)); % pt L - luminozitate
        imA = uint8(im(:, :, 2) + 128); % pt A - axe de culoare
        imB = uint8(im(:, :, 3) + 128); % pt B - axe de culoare

        imLabModif = cat(3, imL, imA, imB);

        [~ ,fName, fileExt] = fileparts(oldFileName);
        newFileNameL = strcat(pathLab, 'L\', fName, fileExt); % path cu imaginele convertite in lab pe foaie L
        newFileNameA = strcat(pathLab, 'A\', fName, fileExt); % path cu imaginele convertite in lab pe foaie A
        newFileNameB = strcat(pathLab, 'B\', fName, fileExt); % path cu imaginele convertite in lab pe foaie B
        newFileNameLabModif = strcat(pathLab, 'Lab128\', fName, fileExt); % path cu imaginele convertite in lab modificat
        imwrite(imL, newFileNameL);
        imwrite(imA, newFileNameA);
        imwrite(imB, newFileNameB);
        imwrite(imLabModif, newFileNameLabModif);
    end
end

function conversionFromRgbToHsv(imdsRgb, pathHsv)
    % Functie utilizata la conversia imaginilor din rgb in hsv. Canalele color space-ul HSV au valori in intervalul
    % (0, 1) asadar este necesara inmultirea fiecarui canal H, S si V cu 128 inainte 
    % de a le stoca local pt a evita pierderea pixelilor. 
    for i = 1 : numel(imdsRgb.Files)
        oldFileName = imdsRgb.Files{i};
           
        im = imread(oldFileName);
        im = rgb2hsv(im);
        
        imH = uint8(im(:, :, 1)*255); % pt H
        imS = uint8(im(:, :, 2)*255); % pt S
        imV = uint8(im(:, :, 3)*255); % pt V
        
        imHSVModif = cat(3, imH, imS, imV);

        [~ ,fName, fileExt] = fileparts(oldFileName);
        newFileNameH = strcat(pathHsv, 'H\', fName, fileExt); % path cu imaginele convertite in hsv pe foaie h
        newFileNameS = strcat(pathHsv, 'S\', fName, fileExt); % path cu imaginele convertite in hsv pe foaie s
        newFileNameV = strcat(pathHsv, 'V\', fName, fileExt); % path cu imaginele convertite in hsv pe foaie v
        newFileNameHsvModif = strcat(pathHsv, 'HsvMod\', fName, fileExt); % path cu imaginele convertite in hsv modificat
        imwrite(imH, newFileNameH);
        imwrite(imS, newFileNameS);
        imwrite(imV, newFileNameV);
        imwrite(imHSVModif, newFileNameHsvModif);
    end
end

function storesLocallyImagesGenerated(imdsGrey, imds_target, FileName, netTrained)
    % functie utilizata pentru stocarea locala a imaginilor generate. 
    for i = 1 : numel(imdsGrey.Files)
        % obtin imaginea generata prin activation si o pun intr-un array
        outputFromNetwork = activations(netTrained, imread(imdsGrey.Files{i}),"ssimL1Loss");
        
        [~, fName, fileExt] = fileparts(imdsGrey.Files{i});
        newFileName = strcat(imds_target, "\" + FileName + "\", fName, fileExt);
        imwrite(uint8(outputFromNetwork), newFileName);
    end
end

function storesLocallyImagesGenerated_ch(imdsGrey, imds_target, FileName, netTrained_ch1, netTrained_ch2, netTrained_ch3, typeFormat)
    % functie utilizata pentru a obtine path-ul cu imagini generate dupa antrenarea retelei
    % SALVEAZA LOCAL O IMAGINE CONCATENATA DE TIP RGB
    for i = 1 : numel(imdsGrey.Files)
        % generare imagini
        outputFromNetwork_ch1 = activations(netTrained_ch1, imread(imdsGrey.Files{i}), "ssimL1Loss");
        outputFromNetwork_ch2 = activations(netTrained_ch2, imread(imdsGrey.Files{i}), "ssimL1Loss");
        outputFromNetwork_ch3 = activations(netTrained_ch3, imread(imdsGrey.Files{i}), "ssimL1Loss");

        if typeFormat == "Lab"
            % scad valorile pt a ajunge in forma dinaintea stocari locale
            outputFromNetwork_ch2_New = outputFromNetwork_ch2 - 128;
            outputFromNetwork_ch3_New = outputFromNetwork_ch3 - 128;
            % salvez local tot o imagine RGB formata prin concatenarea celor 3 canale Lab si convertita in rgb
            outputFromNetwork = cat(3, outputFromNetwork_ch1, outputFromNetwork_ch2_New, outputFromNetwork_ch3_New);
            outputFromNetworkRgb = lab2rgb(outputFromNetwork);
        elseif typeFormat == "HSV"
            % impart valorile pt a ajunge in forma dinaintea stocari locale
            outputFromNetwork_ch1_New = outputFromNetwork_ch1/255;
            outputFromNetwork_ch2_New = outputFromNetwork_ch2/255;
            outputFromNetwork_ch3_New = outputFromNetwork_ch3/255;
             % salvez local tot o imagine RGB formata prin concatenarea celor 3 canale hsv si convertita in rgb
            outputFromNetwork = cat(3, outputFromNetwork_ch1_New, outputFromNetwork_ch2_New, outputFromNetwork_ch3_New);
            outputFromNetworkRgb = hsv2rgb(outputFromNetwork);
        end

        [~, fName, fileExt] = fileparts(imdsGrey.Files{i});
        newFileName = strcat(imds_target, "\" + FileName + "\", fName, fileExt);
        imwrite(uint8(outputFromNetworkRgb*255), newFileName);
    end
end

function targetFormation_residues(imdsOrg, imdsGen, imds_target, FileName, typeFormat)
    % imdsOrg este setul de date tinta, imdsGen este setul de img generate, imds_target este un cale unde pun noile img, \
    % FileName este numele setului de date, typeFormat este tipul formatului de culori
    for i = 1 : numel(imdsOrg.Files)
        if typeFormat == "HSV" || typeFormat == "RGB"
            % incarc imaginile in format hsv stocate local => au valori in intervalul (0, 255) ca au fost inmultite cu 255
            im = single(imread(imdsOrg.Files{i}));
            % incarc imaginile generate hsv stocate local => au valori in intervalul (0, 255) ca au fost inmultite cu 255
            imGen = single(imread(imdsGen.Files{i}));
            imRez = (im - imGen + 255)/2;
        elseif typeFormat == "Lab"
            % incarc imaginile in format lab stocate local => au valori pt L in intervalul (0, 100) si pt a* si b* in intervalul (0, 255) ca au fost adunate cu 128
            im = single(imread(imdsOrg.Files{i}));
            % incarc imaginile generate lab stocate local => au valori pt L in intervalul (0, 100) si pt a* si b* in intervalul (0, 255) ca au fost adunate cu 128
            imGen = single(imread(imdsGen.Files{i}));
            imRez = im - imGen;
            imRez(:,:,1) = (imRez(:,:,1) + 100)/2; % => valori intre (0, 100)
            imRez(:,:,2) = (imRez(:,:,2) + 255)/2; % => valori intre (0, 255)
            imRez(:,:,3) = (imRez(:,:,3) + 255)/2;
        end
      
        [~ ,fName, fileExt] = fileparts(imdsOrg.Files{i});
        newFile = strcat(imds_target, "\" + FileName + "\", fName, fileExt);
        imwrite(uint8(imRez), newFile);
    end
end

function targetFormation_residues_ch(imdsOrg, imdsGen, imds_target, FileName, FileName1,  FileName2, FileName3, typeFormat)
    % imdsOrg - contine imaginile nemodificate considerate seturi tinta, imdsGen - este o foaie obtinuta dupa antrenare(foi precum L/A/B/H/S/V), 
    % imds_target - path-ul cu care creez calea unde vreau sa stochez imaginile, FileName - contine numele la noul folder, 
    % typeFormat -  este foaia la care vreau sa-i calculez reziduurile
    for i = 1 : numel(imdsOrg.Files)
        im = imread(imdsOrg.Files{i});
        if typeFormat == "Lab"
            % incarc imaginile in format lab stocate local => au valori pt L in intervalul (0, 100) si pt a* si b* in intervalul (0, 255) ca au fost adunate cu 128
            im = single(imread(imdsNemod.Files{i}));
            % incarc imaginile generate lab stocate local => au valori pt L in intervalul (0, 100) si pt a* si b* in intervalul (0, 255) ca au fost adunate cu 128
            imGen = single(imread(imdsGen.Files{i}));
            imRez = (im - imGen + 255)/2;
            imRez(:,:,1) = (imRez(:,:,1) + 100)/2; % => valori intre (0, 100)
            imRez(:,:,2) = (imRez(:,:,2) + 255)/2; % => valori intre (0, 255)
            imRez(:,:,3) = (imRez(:,:,3) + 255)/2;
        elseif typeFormat == "HSV" || "RGB"
            % incarc imaginile in format hsv stocate local => au valori in intervalul (0, 255) ca au fost inmultite cu 255
            im = single(imread(imdsOrg.Files{i}));
            % incarc imaginile generate hsv stocate local => au valori in intervalul (0, 255) ca au fost inmultite cu 255
            imGen = single(imread(imdsGen.Files{i}));
            imRez = (im - imGen + 255)/2;
        end

        % Stochez local reziduurile gasite
        [~ ,fName, fileExt] = fileparts(imdsRgb.Files{i});
        newFileH = strcat(imds_target, "\" + FileName + "\", fName, fileExt);
        imwrite(im, newFileH);
        newFileName = strcat(imds_target, "\" + FileName + "\", fName, fileExt);
        imwrite(d, newFileName);
        newFile1 = strcat(imds_target, "\" + FileName1 + "\", fName, fileExt);
        imwrite(uint8(imRez(:, :, 1)), newFile1);
        newFile2 = strcat(imds_target, "\" + FileName2 + "\", fName, fileExt);
        imwrite(uint8(imRez(:, :, 2)), newFile2);
        newFile3 = strcat(imds_target, "\" + FileName3 + "\", fName, fileExt);
        imwrite(uint8(imRez(:, :, 3)), newFile3);
    end
end

