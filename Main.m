%%%%% Main Function of Code %%%%%

clc;
clear all;
close all;

%% Getting the image data %%
cur_dir = pwd;
choice = 0;
filepath = 'Task01_BrainTumour\Task01_BrainTumour\imagesTr\';
FourD = 'y';

while choice ~= 1 && choice ~= 2 && choice ~= 3

    fprintf('Welcome to the Brain Tumor Segmentation Project! :D\n')
    choice = input('This program runs from the Medical Segmentation Decathlon Brain Tumor Data Set!\n\nDo you want to:\n1) Use a random scan\n2) Pick a scan from the dataset?\n3) Use your own .nii.gz Brain Scan! (May not work well)\n');
    if choice == 1
        scan_num = randi([1, 484]);
        fprintf('You got MRI scan #%d!', scan_num)
        if scan_num < 10
            scan_num = strcat('00', int2str(scan_num));
        elseif scan_num < 100
            scan_num = strcat('0', int2str(scan_num));
        else
            scan_num = int2str(scan_num);
        end
        filepath = strcat(filepath,'BRATS_', scan_num,'.nii.gz');
    elseif choice == 2
        scan_num = input('What scan number do you wanna use? (Between 1 and 484) ');
        if scan_num < 10
            scan_num = strcat('00', int2str(scan_num));
        elseif scan_num < 100
            scan_num = strcat('0', int2str(scan_num));
        else
            scan_num = int2str(scan_num);
        end
        filepath = strcat(filepath,'BRATS_', scan_num,'.nii.gz');
    elseif choice == 3
        fprintf('(Don''t include the name of your .nii.gz file) Example: Task01_BrainTumour\\Task01_BrainTumour\\imagesTr\\')
        filepath = input('\nPlease enter the filepath to your .nii.gz Brain MRI Scan!\n', 's');
        fprintf('Example: BRATS_222.nii.gz\n')
        filename = input('What is the name of your .nii.gz file?\n ', 's');
        FourD = input('Is your .nii.gz file 4D? (y/n) ', 's');
        filepath = strcat(filepath, filename);
        fprintf('Your .nii.gz file is located at: %s\n', filepath);
        fprintf('\n***For this project, the front of the brain should be towards the right of the image.***\n')
    else
        fprintf('Not a valid choice :(\nPlease choose a valid number.')
    end
end


%% Turn from .nii.gz to .nii %%

%'M:\Comp Apps\Final_Project\Task01_BrainTumour\'

if choice == 1 || choice == 2
    nii_Brain = gunzip(filepath, 'Task01_BrainTumour\');
    fprintf('\n.nii file stored in the folder "Task01_BrainTumour\"!\n')
    path = 'Task01_BrainTumour\';
    file = filepath(end-15:end-3);
    fileroot = filepath(end-15:end-7);
else
    nii_Brain = gunzip(filepath, pwd);
    fprintf('\n.nii file stored in the present working directory!\n')
    path = pwd;
    file = filename(1:end-3);
    fileroot = filename(1:end-7);
end
    

%% Turn from .nii to .png %%

nii2png_BT(path, file, choice)

%% Poor but Fast Segmentation to choose which slice to use %%
cd(cur_dir);

niftiscan = niftiread(filepath);

if FourD == 'y'
    [~, ~, num_slices, dimes] = size(niftiscan);
    slice_name = SliceDecider(fileroot, num_slices, FourD, path);
else
    [~, ~, num_slices] = size(niftiscan);
    slice_name = SliceDecider(fileroot, num_slices, FourD, path);
end

figure
imshow(slice_name)
title('MRI Brain Scan Slice')

%% Watershed Segmentation on 'Best' Slice %%

I = imread(slice_name);
I = im2gray(I);

cd(cur_dir);

gmag = imgradient(I);
% imshow(gmag,[])
% title("Gradient Magnitude")

L = watershed(gmag);
Lrgb = label2rgb(L);
% imshow(Lrgb)
% title("Watershed Transform of Gradient Magnitude")

se = strel("disk",5);
Io = imopen(I,se);
% imshow(Io)
% title("Opening")

Ie = imerode(I,se);
Iobr = imreconstruct(Ie,I);
% imshow(Iobr)
% title("Opening-by-Reconstruction")

Ioc = imclose(Io,se);
% imshow(Ioc)
% title("Opening-Closing")

Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
% imshow(Iobrcbr)
% title("Opening-Closing by Reconstruction")

fgm = imregionalmax(Iobrcbr);
% imshow(fgm)
% title("Regional Maxima of Opening-Closing by Reconstruction")

I2 = labeloverlay(I,fgm);
% imshow(I2)
% title("Regional Maxima Superimposed on Original Image")

se2 = strel(ones(5,5));
fgm2 = imclose(fgm,se2);
fgm3 = imerode(fgm2,se2);

fgm4 = bwareaopen(fgm3,20);
I3 = labeloverlay(I,fgm4);
% imshow(I3)
% title("Modified Regional Maxima Superimposed on Original Image")

bw = imbinarize(Iobrcbr);
% imshow(bw)
% title("Thresholded Opening-Closing by Reconstruction")

D = bwdist(bw);
DL = watershed(D);
bgm = DL == 0;
% imshow(bgm)
% title("Watershed Ridge Lines")

gmag2 = imimposemin(gmag, bgm | fgm4);

L = watershed(gmag2);

labels = imdilate(L==0,ones(3,3)) + 2*bgm + 3*fgm4;
% I4 = labeloverlay(I,labels);
% imshow(I4)
% title("Markers and Object Boundaries Superimposed on Original Image")

Lrgb = label2rgb(L,"jet","w","shuffle");
figure(2)
imshow(Lrgb)
title("Colored Watershed Label Matrix")

figure(3)
imshow(I)
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.3;
title("Colored Labels Superimposed Transparently on Original Image")

%% Binary Mask of Entire Brain %%

BW_I = imbinarize(I, 0.01);
BW_Edges = edge(BW_I, "sobel");
seT = strel("disk",1);
BW_Edges2 = imdilate(BW_Edges,seT);
L(L==0) = L(1,1);
num_seg = unique(L);

[w, l] = size(BW_I);

%% Choose the Segment that has the Tumor %%

tumor_check = zeros(1, length(num_seg));

disp('Processing, this may take a several minutes...')


skip = 0;
for ii = min(num_seg):1:max(num_seg)
    skip = 0;
    for jj = 1:1:w
        if skip == 0
            for kk = 1:1:l
                if BW_Edges2(jj,kk) == 1 && L(jj,kk) == ii 
                    skip = 1;
                    tumor_check(1,ii) = 0;
                    break
                elseif L(jj,kk) == L(1,1) && L(jj,kk) == ii
                    skip = 1;
                    tumor_check(1,ii) = 0;
                    break
                else
                    tumor_check(1,ii) = ii;
                end
            end
        else
            break
        end
    end
end


disp('Done!')

%Tumor Check: 1 = Label has Tumor, 0 = Label Doesn''t Have Tumor'
disp(tumor_check)

if sum(tumor_check) > 0    %ismember(1, tumor_check)
    fprintf('Tumor Detected!!!\n')
else
    fprintf('No Tumor :(\n')
end

%% Binary Image of Tumor %%
Tumor_Bin = zeros(w, l);

for ii = min(num_seg):1:max(num_seg)
    if tumor_check(ii) ~= 0
        Tumor_Bin(L==ii) = 1;
    end
end

%% Alignment / Image Registration %%

%%adjusting the size,center, and orientation of our reference image to be the same
%%as our scans
Reference = imread('Reference_scan.png');
Reference = imrotate(Reference,-90);
testScan = I;
%shows the positions of the two scans relative to eachother
meh = imshow(imfuse(Reference,testScan)); %this is wayyyyy off, first comparison
%adjust the size of the reference
better = imresize(Reference,[w, l]);
imshow(imfuse(testScan,better)); %this is much much better, this line just compares
imwrite(better,"adjusted_ref.png");

% Making the binary of the lobes

ref = imread('COLOR.png');
ref = imresize(ref, [w, l]); %ensure that they're the same size
imwrite(ref,'Color_actual.png');

%%red = frontal lobe%%
[r,c,color] = size(ref) ;
binaryFrontal = false(r,c) ;
for i = 1:1:r
    for j = 1:1:c
        R = ref(i,j,1);
        G = ref(i,j,2);
        B = ref(i,j,3);
        if  R == 255 && B == 0 && G == 0
            binaryFrontal(i,j) = 1;
        end
    end
end

%imshow(binaryFrontal) %perfect :3
%%green = temporal lobe%%
binaryParietal = false(r,c) ;
for i = 1:1:r
    for j = 1:1:c
        R = ref(i,j,1);
        G = ref(i,j,2);
        B = ref(i,j,3);
        if  R == 0 && B == 0 && G == 255
            binaryParietal(i,j) = 1;
        end
    end
end

%imshow(binaryTemporal)
%%blue = occipital lobe%%
binaryOccipital = false(r,c) ;
for i = 1:1:r
    for j = 1:1:c
        R = ref(i,j,1);
        G = ref(i,j,2);
        B = ref(i,j,3);
        if  R == 0 && B == 255 && G == 0
            binaryOccipital(i,j) = 1;
        end
    end
end

%imshow(binaryOccipital)
%%white = temporal lobe%%
binaryTemporal= false(r,c) ;
for i = 1:1:r
    for j = 1:1:c
        R = ref(i,j,1);
        G = ref(i,j,2);
        B = ref(i,j,3);
        if  R == 255 && B == 255 && G == 255
            binaryTemporal(i,j) = 1;
        end
    end
end

%imshow(binaryParietal)
imwrite(binaryFrontal,'Frontal.png');
imwrite(binaryParietal,'Parietal.png');
imwrite(binaryOccipital,'Occipital.png');
imwrite(binaryTemporal,'Temporal.png');

Frontalcount = 0;
Parietalcount = 0;
Occipitalcount = 0;
Temporalcount = 0;

%% Calculating Where The Tumor Is %%

TumorArea = bwarea(Tumor_Bin);

for ii = 1:1:w
    for jj = 1:1:l
        if Tumor_Bin(ii,jj) == 1 && binaryFrontal(ii,jj) == 1
            Frontalcount = Frontalcount + 1;
        elseif Tumor_Bin(ii,jj) == 1 && binaryParietal(ii,jj) == 1
            Parietalcount = Parietalcount + 1;
        elseif Tumor_Bin(ii,jj) == 1 && binaryOccipital(ii,jj) == 1
            Occipitalcount = Occipitalcount + 1;
        elseif Tumor_Bin(ii,jj) == 1 && binaryTemporal(ii,jj) == 1
            Temporalcount = Temporalcount + 1;
        end
    end
end

Fper = (Frontalcount/TumorArea)*100;
Pper = (Parietalcount/TumorArea)*100;
Oper = (Occipitalcount/TumorArea)*100;
Tper = (Temporalcount/TumorArea)*100;

Per = [Fper, Pper, Oper, Tper];

Mper = max(Per);

if Mper == Fper
    fprintf('%.2f%% of the tumor is in the frontal lobe!', Mper)
elseif Mper == Pper
    fprintf('%.2f%% of the tumor is in the parietal lobe!', Mper)
elseif Mper == Oper
    fprintf('%.2f%% of the tumor is in the occipital lobe!', Mper)
elseif Mper == Tper
    fprintf('%.2f%% of the tumor is in the temporal lobe!', Mper)
end


figure(4)
imshow(imfuse(Tumor_Bin,ref))