clear all
close all
clc

%% parametri
[dataDir, ~, ~] = fileparts(mfilename('fullpath')); % incarca in dataDir calea curenta a script-ului
pathTestImIn = datastore(strcat(dataDir, "\test\testSize\in")); % imagini gri de testare
pathTestImOut_org = datastore(strcat(dataDir, "\test\testSize\out")); % imagini color de testare
afisare_toate_img_generate = false;
testare_grey_to_rgb = false;
testare_HSV = false;
testare_HSV_pe_foi = false;
testare_LAB = false;
testare_LAB_pe_foi = false;
test_res_grey_to_rgb = false;
test_res_grey_to_lab = false;
test_res_grey_to_lab_pe_foi = false;
test_res_grey_to_hsv = true;
nr = numel(pathTestImIn.Files);
x = randi(nr);

if afisare_toate_img_generate == true
    afisare_rezultate(pathTestImIn, pathTestImOut_org);
end

%% testare imagini din tonuri de gri in color
if testare_grey_to_rgb == true
    % incarca retelele antrenate
    % load('GreyToRgb_100epoci_antrenareRetea-2023-03-10-18-40-27.mat'); % 100 epoci, 20mB, retea preantrenata
    load('greyToRgb_100+50_epoci_CuValFreq250--2023-04-17-21-49-03.mat'); % 150 epoci, 20mB
    netOldAntre = netTrained;

    outputFromNetwork_rgb = [];
    for i = 1 : nr
        % obtin imaginea generata prin activation si o pun intr-un array
        outputFromNetwork_rgb = cat(4,outputFromNetwork_rgb, activations(netOldAntre,imread(pathTestImIn.Files{i}),"ssimL1Loss"));
    end

    % calculez calitatea img generate cu img color originale
    [qErrPatratica_RGB, qSSIM_RGB, qNIQE_RGB, qPSNR_RGB] = imagesQuality(pathTestImOut_org, uint8(outputFromNetwork_rgb));
    qErrPFinal_RGB = (qErrPatratica_RGB(:, :, 1) + qErrPatratica_RGB(:, :, 2) + qErrPatratica_RGB(:, :, 3)) / 3;

    subplot(1, 2, 1);
    imshow(imread(pathTestImOut_org.Files{x})); title('img originala');
    subplot(1, 2, 2);
    imshow(uint8(outputFromNetwork_rgb(:, :, :, x))); title('img generata GREY-RGB');
end

%% testare LAB
if testare_LAB == true
    % incarca imaginile pentru testare convertite in lab 
    imds_LAB = imageDatastore(strcat(dataDir,"\test\testLab\outLab128"));
    % load('antrenareReteaGreyToLAB_100Epoci_FolosindLABOut128_CuValFreq250-2023-03-27-01-14-01.mat'); % 100 de epoci
    load('antrenareReteaGreyToLAB_100_50Epoci_FolosindLABOut128_CuValFreq250-2023-04-14-17-43-00.mat'); % 150 de epoci
    netOld_LAB_Mod = netTrained;

    outputFromNetwork_LAB_Mod = [];
    for i = 1 : nr
        % obtin imaginea generata prin activation si o pun intr-un array
        outputFromNetwork_LAB_Mod = cat(4,outputFromNetwork_LAB_Mod, activations(netOld_LAB_Mod,imread(pathTestImIn.Files{i}),"ssimL1Loss"));
        % pt a o putea compara cu img originala scad cu 128 foile A si B din array
        outputFromNetwork_LAB_Mod(:, :, 2, i) = outputFromNetwork_LAB_Mod(:, :, 2, i) - 128;
        outputFromNetwork_LAB_Mod(:, :, 3, i) = outputFromNetwork_LAB_Mod(:, :, 3, i) - 128;
    end
    % convertesc in rgb, dar inmultite cu 255 ca valorile dupa conversie sa fie intre 0 si 255 (inainte erau valori intre 0 si 1)
    outputFromNetwork_LAB_Mod1  = lab2rgb(outputFromNetwork_LAB_Mod)*255;

    % calculez calitatea img generate cu img color originale
    [qErrPatratica_LAB, qSSIM_LAB, qNIQE_LAB, qPSNR_LAB] = imagesQuality(pathTestImOut_org, uint8(outputFromNetwork_LAB_Mod1));
    qErrPFinal_LAB = (qErrPatratica_LAB(:, :, 1) + qErrPatratica_LAB(:, :, 2) + qErrPatratica_LAB(:, :, 3)) / 3;

    subplot(1, 2, 1);
    imshow(imread(pathTestImOut_org.Files{x})); title('img originala');
    subplot(1, 2, 2);
    imshow(uint8(outputFromNetwork_LAB_Mod1(:, :, :, x))); title('img generata GREY-LAB');
end

%% LAB pe foi
if testare_LAB_pe_foi == true
    % FOAIE L
    imdsL = imageDatastore(strcat(dataDir, "\test\testLab\outL")); % contine imagini lab - fereastara L
    load('antrenareFoaie_L_100Epoci-2023-03-24-17-17-52.mat'); % pe foaia L, 100 epoci, FARA VALIDARE
    netOldL100 = netTrained;
    load('antrenareFoaie_L_100_50Epoci-CuValFreq250-2023-04-15-11-32-20.mat'); % 150 de epoci
    netOldL100_50 = netTrained;
    % FOAIA A
    imdsA = imageDatastore(strcat(dataDir, "\test\testLab\outA")); % contine imagini lab - fereastara A
    load('GreyToLab_Foaia_A_100epoci_net_checkpoint__50000__2023_03_25__05_22_48.mat'); % 100 de epoci
    netOldA100 = net;
    load('antrenareFoaie_A_150EpociCuValFreq250-2023-03-25-09-43-18.mat'); % 150 epoci, foaia A, fara sa scad cu 128
    netOldA150 = netTrained;
    load('antrenareFoaie_A_150+150EpociCuValFreq250-2023-03-26-09-34-23.mat'); % foaia A, 300 epoci
    netOldA150_150=netTrained;

    % FOAIA B
    imdsB = imageDatastore(strcat(dataDir,"\test\testLab\outB")); % contine imagini lab - fereastara B
    load('greyToLab_foaia_B_100epoci_net_checkpoint__50000__2023_03_25__02_33_48.mat');
    netOldB100 = net;
    load('net_checkpoint__75000__2023_03_25__09_52_57-foaiaB-150Ep.mat'); % 150 de epoci, foaia B checkpoint
    netOldB150_ch = net;
    load('antrenareFoaie_B_150+150EpociCuValFreq250-2023-03-25-23-48-15.mat'); % 300 epoci, foaia B
    netOldB150_150 = netTrained;

    outputFromNetworkL = []; outputFromNetworkA = []; outputFromNetworkB = []; outputFromNetworkLAB = []; labToRgb = []; outputFromNetworkLABToRgb = [];
    for i = 1 : nr
        % pun in array-uri imaginile generate de retea pe fiecare foaie
        outputFromNetworkL  = cat(3,outputFromNetworkL,activations(netOldL100_50,imread(pathTestImIn.Files{i}),"ssimL1Loss"));
        outputFromNetworkA  = cat(3,outputFromNetworkA,activations(netOldA150,imread(pathTestImIn.Files{i}),"ssimL1Loss"));
        outputFromNetworkB  = cat(3,outputFromNetworkB,activations(netOldB150_ch,imread(pathTestImIn.Files{i}),"ssimL1Loss"));
    end

    outputFromNetworkLABToRgb  = cat(4,outputFromNetworkLABToRgb, outputFromNetworkL, outputFromNetworkA-128, outputFromNetworkB-128);
    outputFromNetworkLABToRgb = permute(outputFromNetworkLABToRgb, [1 2 4 3]);
    labToRgb = cat(4,labToRgb,lab2rgb(outputFromNetworkLABToRgb)*255); % *255 pentru ca fiecare pixel sa aiba valori intre 0 si 255, nu intre 0 si 1 cum returneaza lab2rgb
    
    
    [qErrPatratica_L, qSSIM_L, qNIQE_L, qPSNR_L] = imagesQuality(imdsL, uint8(outputFromNetworkL));
    [qErrPatratica_A, qSSIM_A, qNIQE_A, qPSNR_A] = imagesQuality(imdsA, uint8(outputFromNetworkA));
    [qErrPatratica_B, qSSIM_B, qNIQE_B, qPSNR_B] = imagesQuality(imdsB, uint8(outputFromNetworkB));
    [qErrPatratica_LAB_rgb, qSSIM_LAB_rgb, qNIQE_LAB_rgb, qPSNR_LAB_rgb] = imagesQuality(pathTestImOut_org, uint8(labToRgb));
    qErrPatratica_LAB_rgb = (qErrPatratica_LAB_rgb(:, :, 1) + qErrPatratica_LAB_rgb(:, :, 2) + qErrPatratica_LAB_rgb(:, :, 3))/3;
    
    subplot(1, 2, 1);
    imshow(imread(pathTestImOut_org.Files{x})); title('img originala');
    subplot(1, 2, 2);
    imshow(uint8(labToRgb(:, :, :, x))); title('img generata');

end

%% testare HSV
if testare_HSV == true
    load('antrenareHSV_150Epoci_FolosindHSVModificat-2023-03-30-22-59-44.mat'); % 150 de epoci
    netOld_HSV_Mod = netTrained;

    outputFromNetwork_HSV_Mod = []; hsvToRgb = [];
    for i = 1 : nr
        % obtin imaginea generata prin activation si o pun intr-un array
        outputFromNetwork_HSV_Mod = cat(4,outputFromNetwork_HSV_Mod, activations(netOld_HSV_Mod,imread(pathTestImIn.Files{i}),"ssimL1Loss"));
        % pt a o putea compara cu img originala inmultesc cu 255 toate foile
        outputFromNetwork_HSV_Mod(:, :, :, i) = outputFromNetwork_HSV_Mod(:, :, :, i)/255;
        % dupa care convertesc in rgb, dar inmultite cu 255 ca valorile dupa conversie sa fie intre 0 si 255 (inainte erau valori intre 0 si 1)
        hsvToRgb = cat(4, hsvToRgb, hsv2rgb(outputFromNetwork_HSV_Mod(:, :, :, i)) * 255);
    end
    % calculez calitatea img generate convertite in rgb cu img color originale
    [qErrPatratica_hsv, qSSIM_hsv, qNIQE_hsv, qPSNR_hsv] = imagesQuality(pathTestImOut_org, uint8(hsvToRgb));
    qErrPFinal_hsv = (qErrPatratica_hsv(:, :, 1) + qErrPatratica_hsv(:, :, 2) + qErrPatratica_hsv(:, :, 3)) / 3;

    subplot(1, 2, 1);
    imshow(imread(pathTestImOut_org.Files{x})); title('img originala');
    subplot(1, 2, 2);
    imshow(uint8(hsvToRgb(:, :, :, x))); title('img generata GREY-HSV');

end
%% HSV pe canale
if testare_HSV_pe_foi == true
    % FOAIE H
    imdsH = imageDatastore(strcat(dataDir, "\test\testHsv\outH")); % imagini hsv stocate local - fereastara h
    load('antrenareFoaie_H_150 epoci-2023-03-29-14-30-30.mat'); % foaia h, 150 epoci
    netOldH150 = netTrained;
    load('antrenareFoaie_H_150+150epoci-2023-04-04-23-13-04.mat');
    netOldH150_150 = netTrained;
    
    % FOAIE S
    imdsS = imageDatastore(strcat(dataDir, "\test\testHsv\outS")); % imagini hsv stocate local - fereastara S
    load('antrenareFoaie_S_150epoca-2023-03-29-05-53-38.mat'); % foaia S, 150 epoci
    netOldS150 = netTrained;
    load('antrenareFoaie_S_150+150epoci_net_checkpoint__75000__2023_04_05__02_17_02.mat');
    netOldS150_150 = net;
    
    % FOAIE V
    imdsV = imageDatastore(strcat(dataDir,"\test\testHsv\outV")); % imagini hsv stocate local - fereastara V
    load('antrenareFoaie_V_150epoca-2023-03-29-17-02-32.mat'); % 150 epoci, foaia V
    netOldV150 = netTrained;
    
    outputFromNetworkH = []; outputFromNetworkS = []; outputFromNetworkV = []; outputFromNetwork_HSV = []; hsvToRgb = [];
    for i = 1 : nr
        % pun in array-uri diferite imaginile generate de fiecare retea
        outputFromNetworkH  = cat(3,outputFromNetworkH,activations(netOldH150,imread(pathTestImIn.Files{i}),"ssimL1Loss"));
        outputFromNetworkS  = cat(3,outputFromNetworkS,activations(netOldS150,imread(pathTestImIn.Files{i}),"ssimL1Loss"));
        outputFromNetworkV  = cat(3,outputFromNetworkV,activations(netOldV150,imread(pathTestImIn.Files{i}),"ssimL1Loss"));
    end
    % concatenez array-urile de imagini generate fara a imparti la 255 sa pot compara cu imds-urile HSV_modificate local
    outputFromNetwork_HSV  = cat(4,outputFromNetwork_HSV, outputFromNetworkH, outputFromNetworkS, outputFromNetworkV);
    % concatenez array-urile de imagini generate, dar impartite la 255 pt a converti imaginile la rgb
    outputFromNetwork_HSV_for_RGB = outputFromNetwork_HSV/255;
    % pemutare astfel incat sa aiba forma linie*coloana*canale*nrImg 
    outputFromNetwork_HSV_for_RGB = permute(outputFromNetwork_HSV_for_RGB, [1 2 4 3]);

    for i = 1 : nr
        % le convertesc in rgb si inmultesc cu 255 sa fie sa fie in acelasi interval cu imaginile rgb originale
        hsvToRgb = cat(4,hsvToRgb,(hsv2rgb(outputFromNetwork_HSV_for_RGB(:, :, :, i)))*255);
    end

    % apelez functia de calculare a calitatii pe fiecare foaie, pe hsv si rgb
    [qErrPatratica_H, qSSIM_H, qNIQE_H, qPSNR_H] = imagesQuality(imdsH, uint8(outputFromNetworkH));
    [qErrPatratica_S, qSSIM_S, qNIQE_S, qPSNR_S] = imagesQuality(imdsS, uint8(outputFromNetworkS));
    [qErrPatratica_V, qSSIM_V, qNIQE_V, qPSNR_V] = imagesQuality(imdsV, uint8(outputFromNetworkV));
    [qErrPatratica_HSV_rgb, qSSIM_HSV_rgb, qNIQE_HSV_rgb, qPSNR_HSV_rgb] = imagesQuality(pathTestImOut_org, uint8(hsvToRgb));

    % calculez eroarea patratica totala pe hsv si rgb
    qErrPatratica_HSV_rgb = (qErrPatratica_HSV_rgb(:, :, 1) + qErrPatratica_HSV_rgb(:, :, 2) + qErrPatratica_HSV_rgb(:, :, 3))/3;
    
    subplot(1, 2, 1);
    imshow(imread(pathTestImOut_org.Files{x})); title('img originala');
    subplot(1, 2, 2);
    imshow(uint8(hsv2rgb(outputFromNetwork_HSV_for_RGB(:, :, :, x))*255)); title('img generata');
end

%% adunare reziduuri la imaginile generate din grey to rgb
if test_res_grey_to_rgb == true
    pathTestImOut_Gen = datastore(strcat(dataDir,"\test\Img_Generate_grey_to_rgb")); % img rgb generate
    load('GreyToReziduuri_PlecandDeLa_GreyToRgb_50_epoci_varianta profa-2023-04-07-19-39-34.mat');
    netOld_reziduuri_30 = netTrained;

    outputFromNetwork_reziduuri_30 = []; outputFromNetwork_rez_rgbGen= [];
    for i = 1 : numel(pathTestImOut_Gen.Files)
        rgbGen = single(imread(pathTestImOut_Gen.Files{i}));
        imGray = imread(pathTestImIn.Files{i});
        outputFromNetwork_reziduuri_30 = cat(4, outputFromNetwork_reziduuri_30, activations(netOld_reziduuri_30, imGray,"ssimL1Loss"));
        % restaurez img => adaug reziduurile generate la imaginea generata anterior acum stocata local
        outputFromNetwork_rez_rgbGen = cat(4, outputFromNetwork_rez_rgbGen, uint8(outputFromNetwork_reziduuri_30(:,:,:,i) * 2 - 255 + rgbGen));
    end
    [qErrPatratica, qSSIM, qNIQE, qPSNR] = imagesQuality(pathTestImOut_org, outputFromNetwork_rez_rgbGen);
    
    subplot(1, 2, 1);
    imshow(imread(pathTestImOut_org.Files{x})); title('IMAGINE ORIGINALA');
    subplot(1, 2, 2);
    imshow(uint8(outputFromNetwork_rez_rgbGen(:, :, :, x))); title('Imagine generata grey to rgb + res');
end

%% adunare reziduuri la imaginile generate din grey to lab
if test_res_grey_to_lab == true
    load('GreyToReziduuri_PlecandDeLa_GreyToLab_pe_foi_100_epoci-2023-04-08-13-35-15.mat');
    netOld_reziduuri_30 = netTrained;
    pathTestImOut_lab_Gen = datastore(strcat(dataDir,"\test\Img_Generate_grey_to_lab"));

    outputFromNetwork_reziduuri_30 = []; outputFromNetwork_rez_labGen = [];
    for i = 1 : numel(pathTestImOut_lab_Gen.Files)
        rgbGen_lab = single(imread(pathTestImOut_lab_Gen.Files{i}));
        imGray = imread(pathTestImIn.Files{i});
        outputFromNetwork_reziduuri_30 = cat(4, outputFromNetwork_reziduuri_30, activations(netOld_reziduuri_30, imGray,"ssimL1Loss"));
        outputFromNetwork_rez_labGen = cat(4, outputFromNetwork_rez_labGen, outputFromNetwork_reziduuri_30(:,:,:,i) * 2 - 255 + rgbGen_lab);
    end
    [qErrPatratica, qSSIM, qNIQE, qPSNR] = imagesQuality(pathTestImOut_org, uint8(outputFromNetwork_rez_labGen));
    subplot(1, 2, 1);
    imshow(imread(pathTestImOut_org.Files{x})); title('img orginala');
    subplot(1, 2, 2);
    imshow(uint8(outputFromNetwork_rez_labGen(:, :, :, x))); title('img generate grey to lab + rez');
end

%% adunare reziduuri la imaginile generate din grey to lab pe foi
if test_res_grey_to_lab_pe_foi == true
    % incarca load-urile cu modelele antrenate spre reziduuri
    load('antrenareFoaie_L_100Epoci-2023-03-24-17-17-52.mat'); % pe foaia L, 100 epoci, FARA VALIDARE
    netOldL100 = netTrained;
    load('GreyToReziduuri_PlecandDeLa_GreyTo_Lab_foaia_A_100_epoci-2023-04-13-21-24-29.mat'); % pe foaia A, 100 epoci,
    netOld_reziduuri_foaia_A = netTrained;
    load('GreyToReziduuri_PlecandDeLa_GreyTo_Lab_foaia_B_100_epoci-2023-04-13-21-39-28.mat'); % pe foaia B, 100 epoci,
    netOld_reziduuri_foaia_B = netTrained;

    % datastore-uri cu canalele generate in etapa anterioara
    pathTestImOut_Gen_lab_ch_A = datastore(strcat(dataDir,"\test\Img_Generate_grey_to_lab_pe_foaia_A"));
    pathTestImOut_Gen_lab_ch_B = datastore(strcat(dataDir,"\test\Img_Generate_grey_to_lab_pe_foaia_B"));
    
    arrRez_A = []; arrRez_B = []; arr_L = []; labToRgb = []; 
    for i = 1 : numel(pathTestImOut_Gen_lab_ch_A.Files)
        chA_Gen = single(imread(pathTestImOut_Gen_lab_ch_A.Files{i}))-128; % => interval (-128, 127)
        chB_Gen = single(imread(pathTestImOut_Gen_lab_ch_B.Files{i}))-128;
        imGray = imread(pathTestImIn.Files{i});

        outputFromNetworkL = activations(netOldL100, imGray,"ssimL1Loss"); % se genereaza canalul L
        outputFromNetwork_reziduuri_foaia_A = single(activations(netOld_reziduuri_foaia_A, imGray,"ssimL1Loss")); % => interval (0, 255)
        outputFromNetwork_reziduuri_foaia_B = single(activations(netOld_reziduuri_foaia_B, imGray,"ssimL1Loss"));
        
        arr_L = cat(3, arr_L, outputFromNetworkL);
        arrRez_A = cat(3, arrRez_A, outputFromNetwork_reziduuri_foaia_A);
        arrRez_B = cat(3, arrRez_B, outputFromNetwork_reziduuri_foaia_B);
        
        ch_a_final1 = (outputFromNetwork_reziduuri_foaia_A*2-255) + chA_Gen;
        ch_b_final1 = (outputFromNetwork_reziduuri_foaia_B*2-255) + chB_Gen;

        outputFromNetwork_LAB_rez2 = cat(4, outputFromNetworkL, ch_a_final1, ch_b_final1);
        outputFromNetwork_LAB_rez2 = permute(outputFromNetwork_LAB_rez2, [1 2 4 3]); 

        labToRgb = cat(4,labToRgb,lab2rgb(outputFromNetwork_LAB_rez2)*255); % *255 pentru ca fiecare pixel sa aiba valori intre 0 si 255, nu intre 0 si 1 cum returneaza lab2rgb
        labToRgb = uint8(labToRgb);
    end

    % se calculeaza calitatea imaginiilor corectate cu reziduuri
    [qErrPatratica, qSSIM, qNIQE, qPSNR] = imagesQuality(pathTestImOut_org, uint8(labToRgb));

    subplot(1, 2, 1);
    imshow(imread(pathTestImOut_org.Files{x})); title('img orginala');
    subplot(1, 2, 2);
    imshow((labToRgb(:, :, :, x))); title('img generate grey to lab pe foi + rez pe foi');
end

%% adunare reziduuri la imaginile generate din grey to hsv pe foi
if test_res_grey_to_hsv == true
    load('GreyToReziduuri_100ep_PlecandDeLa_GreyToHsv(cu_150ep)-2023-06-18-03-30-57.mat');
    netOld_reziduuri_30 = netTrained;

    pathTestImOut_Gen_HSV = datastore(strcat(dataDir, "\test\Img_Generate_grey_to_HSV"));
    
    outputFromNetwork_reziduuri_30 = []; outputFromNetwork_rez_hsvGen = [];
    for i = 1 : numel(pathTestImOut_Gen_HSV.Files)
        rgbGen_hsv = single(imread(pathTestImOut_Gen_HSV.Files{i}));
        imGray = imread(pathTestImIn.Files{i});

        % se genereaza reziduurile
        outputFromNetwork_reziduuri_30 = cat(4, outputFromNetwork_reziduuri_30, activations(netOld_reziduuri_30, imGray,"ssimL1Loss"));
        outputFromNetwork_rez_hsvGen = cat(4, outputFromNetwork_rez_hsvGen, outputFromNetwork_reziduuri_30(:,:,:,i) * 2 - 255 + rgbGen_hsv);
    end
    [qErrPatratica, qSSIM, qNIQE, qPSNR] = imagesQuality(pathTestImOut_org, uint8(outputFromNetwork_rez_hsvGen));
    qErrPFinal = (qErrPatratica(:, :, 1) + qErrPatratica(:, :, 2) + qErrPatratica(:, :, 3)) / 3;

    subplot(1, 2, 1);
    imshow(imread(pathTestImOut_org.Files{x})); title('img orginale');
    subplot(1, 2, 2);
    imshow(uint8(outputFromNetwork_rez_hsvGen(:, :, :, x))); title('img generate grey to hsv + rez');
    
end
%% Support functions
function [qErrPatratica, qSSIM, qNIQE, qPSNR] = imagesQuality(imdsOut, outputFromNetworkImg)
    qErrPatratica = 0; qSSIM = 0; qNIQE = 0; qPSNR = 0;
    nr = numel(imdsOut.Files);
    for i = 1 : nr
        if numel(size(outputFromNetworkImg)) == 3
            % pentru cazul in care outputFromNetwork contine o singura foaie, ex L, H, A (lin*col*nrImg)
            outputFromNetwork = outputFromNetworkImg(:,:,i);
        else
            % pentru cazul in care outputFromNetwork contine toate foile, ex LAB, HSV, RGB (lin*col*ch*nrImg)
            outputFromNetwork = outputFromNetworkImg(:,:,:,i);
        end
        im = imread(imdsOut.Files{i});

        % eroare medie patratica
        dif = (outputFromNetwork - im).^2; % diferenta dintre ce am obtinut si ce voiam sa obtin
        qErrPatratica = qErrPatratica + sum(sum(dif))/(128*128)/nr;

        % calitatea imaginii cu ssim
        imgSSIM = ssim(im, outputFromNetwork);
        qSSIM = qSSIM + imgSSIM/nr;

        % calitatea imaginii cu niqe
        imgNiqe = niqe(uint8(outputFromNetwork));
        qNIQE = qNIQE + imgNiqe/nr;

        % calitatea imaginii cu psnr
        imgPsnr = psnr(im, outputFromNetwork);
        qPSNR = qPSNR + imgPsnr/nr;
    end
end
function afisare_rezultate(pathTestImIn, pathTestImOut_org)
    x = randi(numel(pathTestImIn.Files));
    
    %% Etapa1 - testare simultan pe toate canalele
    figure
    % imagine rgb originala
    subplot(1, 4, 1);
    imshow(imread(pathTestImOut_org.Files{x})); title('IMAGINE ORIGINALA');
    % testare imagini cu retea preantrenata
    %load('GreyToRgb_100epoci_antrenareRetea-2023-03-10-18-40-27.mat'); % 100 epoci, 20mB, retea preantrenata
    load('greyToRgb_100+50_epoci_CuValFreq250--2023-04-17-21-49-03.mat');
    netOldAntre = netTrained;
    outputFromNetwork_grey_to_rgb = activations(netOldAntre,imread(pathTestImIn.Files{x}),"ssimL1Loss");
    subplot(1, 4, 2);
    imshow(uint8(outputFromNetwork_grey_to_rgb)); title('img generata GREY-RGB');

    % testare LAB
    %load('antrenareReteaGreyToLAB_100Epoci_FolosindLABOut128_CuValFreq250-2023-03-27-01-14-01.mat'); %100 epoci
    load('antrenareReteaGreyToLAB_100_50Epoci_FolosindLABOut128_CuValFreq250-2023-04-14-17-43-00.mat');
    netOld_LAB_Mod = netTrained;
    outputFromNetwork_grey_to_lab = activations(netOld_LAB_Mod,imread(pathTestImIn.Files{x}),"ssimL1Loss");
    outputFromNetwork_grey_to_lab(:, :, 2) = outputFromNetwork_grey_to_lab(:, :, 2)-128;
    outputFromNetwork_grey_to_lab(:, :, 3) = outputFromNetwork_grey_to_lab(:, :, 3)-128;
    outputFromNetwork_grey_to_lab = lab2rgb(outputFromNetwork_grey_to_lab);
    subplot(1, 4, 3);
    imshow((outputFromNetwork_grey_to_lab)); title('img generata GREY-LAB');

    %testare hsv
    load('antrenareHSV_150Epoci_FolosindHSVModificat-2023-03-30-22-59-44.mat');
    netOld_HSV_Mod = netTrained;
    outputFromNetwork_grey_to_hsv = activations(netOld_HSV_Mod,imread(pathTestImIn.Files{x}),"ssimL1Loss");
    outputFromNetwork_grey_to_hsv =outputFromNetwork_grey_to_hsv/255;
    outputFromNetwork_grey_to_hsv = hsv2rgb(outputFromNetwork_grey_to_hsv);
    subplot(1, 4, 4);
    imshow(outputFromNetwork_grey_to_hsv); title('img generata GREY-HSV');
    %% Etapa2 - separat pe fiecare canal
    figure
    % imagine rgb originala
    subplot(1, 5, 1);
    imshow(imread(pathTestImOut_org.Files{x})); title('IMAGINE ORIGINALA');
    % testare lab pe foi
    %load('antrenareFoaie_L_100Epoci-2023-03-24-17-17-52.mat');% pe foaia L, 100 epoci, FARA VALIDARE
    load('antrenareFoaie_L_100_50Epoci-CuValFreq250-2023-04-15-11-32-20.mat');
    netOldL100 = netTrained;
    load('antrenareFoaie_A_150EpociCuValFreq250-2023-03-25-09-43-18.mat'); %150 epoci, foaia A, fara sa scad cu 128
    netOldA150 = netTrained;
    load('net_checkpoint__75000__2023_03_25__09_52_57-foaiaB-150Ep.mat'); % 150 de epoci, foaia B checkpoint
    netOldB150_ch = net;
    outputFromNetwork_L = activations(netOldL100,imread(pathTestImIn.Files{x}),"ssimL1Loss");
    outputFromNetwork_A_150 = activations(netOldA150,imread(pathTestImIn.Files{x}),"ssimL1Loss");
    outputFromNetwork_B_150 = activations(netOldB150_ch,imread(pathTestImIn.Files{x}),"ssimL1Loss");
    outputFromNetwork_cat = cat(3, outputFromNetwork_L, outputFromNetwork_A_150-128, outputFromNetwork_B_150-128);
    labToRgb_pe_foi_150 = lab2rgb(outputFromNetwork_cat);
    subplot(1, 5, 2);
    imshow((labToRgb_pe_foi_150)); title('img generata GREY-LAB pe foi 150 epoci');
    % 300 epoci
    load('antrenareFoaie_A_150+150EpociCuValFreq250-2023-03-26-09-34-23.mat'); % foaia A, 150+150 epoci
    netOldA150_150=netTrained;
    load('antrenareFoaie_B_150+150EpociCuValFreq250-2023-03-25-23-48-15.mat'); %150+150 epoci, foaia B
    netOldB150_150 = netTrained;
    outputFromNetwork_A_300 = activations(netOldA150_150,imread(pathTestImIn.Files{x}),"ssimL1Loss");
    outputFromNetwork_B_300 = activations(netOldB150_150,imread(pathTestImIn.Files{x}),"ssimL1Loss");
    outputFromNetwork_cat_300 = cat(3, outputFromNetwork_L, outputFromNetwork_A_300-128, outputFromNetwork_B_300-128);
    labToRgb_pe_foi_300 = lab2rgb(outputFromNetwork_cat_300);
    subplot(1, 5, 3);
    imshow((labToRgb_pe_foi_300)); title('img generata GREY-LAB pe foi 300 epoci');
    
    % testare hsv pe foi
    load('antrenareFoaie_V_150epoca-2023-03-29-17-02-32.mat'); % 150 epoci, foaia V
    netOldV150 = netTrained;
    load('antrenareFoaie_H_150 epoci-2023-03-29-14-30-30.mat'); %foaia h, 150 epoci
    netOldH150 = netTrained;
    load('antrenareFoaie_S_150epoca-2023-03-29-05-53-38.mat'); % foaia S, 150 epoci
    netOldS150 = netTrained;
    outputFromNetwork_H_150 = activations(netOldH150,imread(pathTestImIn.Files{x}),"ssimL1Loss");
    outputFromNetwork_S_150 = activations(netOldS150,imread(pathTestImIn.Files{x}),"ssimL1Loss");
    outputFromNetwork_V_150 = activations(netOldV150,imread(pathTestImIn.Files{x}),"ssimL1Loss");
    outputFromNetwork_cat_HSV_150 = cat(3, outputFromNetwork_H_150/255, outputFromNetwork_S_150/255, outputFromNetwork_V_150/255);
    HSVToRgb_pe_foi_150 = hsv2rgb(outputFromNetwork_cat_HSV_150);
    subplot(1, 5, 4);
    imshow((HSVToRgb_pe_foi_150)); title('img generata GREY-HSV pe foi 150 epoci');

    load('antrenareFoaie_H_150+150epoci-2023-04-04-23-13-04.mat');
    netOldH150_150 = netTrained;
    load('antrenareFoaie_S_150+150epoci_net_checkpoint__75000__2023_04_05__02_17_02.mat');
    netOldS150_150 = net;
    outputFromNetwork_H_300 = activations(netOldH150_150,imread(pathTestImIn.Files{x}),"ssimL1Loss");
    outputFromNetwork_S_300 = activations(netOldS150_150,imread(pathTestImIn.Files{x}),"ssimL1Loss");
    outputFromNetwork_cat_HSV_300 = cat(3, outputFromNetwork_H_300/255, outputFromNetwork_S_300/255, outputFromNetwork_V_150/255);
    HSVToRgb_pe_foi_300 = hsv2rgb(outputFromNetwork_cat_HSV_300);
    subplot(1, 5, 5);
    imshow((HSVToRgb_pe_foi_300)); title('img generata GREY-HSV pe foi 300 epoci');

end