clc;clear;close all;
ageList = [26, 27, 28, 21, 32, 29, 31, 31, 26, 29, ...
        27, 23, 33, 37, 31, 36, 29, 29, 23, 23, ...
        19, 21, 27, 22, 24, 25, 22, 22, 24, 25, ...
        18, 20, 24, 18, 25, 24, 22, 24, 18, 16, ...
        34, 30, 40, 36, 41, 35, 31];
ageMin = min(ageList);
ageMax = max(ageList);
ageMean = mean(ageList);
ageSD = std(ageList);
load('E:\IEDpropagation\step2_fiberTrack\whole_brain_Yeo2011_7Networks_qa.mat');
qa = connectivity;
load('E:\IEDpropagation\step2_fiberTrack\whole_brain_Yeo2011_7Networks_ncount.mat')
ncount = connectivity;
plot(qa(:), log(ncount(:)), '*')
idx = find(ncount>0);
[r,  p] = corr(qa(idx), log(ncount(idx)));
load('whole_brain_Yeo2011_17Networks_MNI152NLin2009cAsym_1mm_LPS_qa.mat')
qa = connectivity;
load('whole_brain_Yeo2011_17Networks_MNI152NLin2009cAsym_1mm_LPS_ncount.mat')
ncount = connectivity;
figure;
subplot(2, 2, 1)
imagesc(log(ncount));
title('ncount');
subplot(2, 2, 2);
imagesc(qa);
title('QA');
subplot(2, 2, 3);
plot(log(ncount(:)), qa(:), '*')
ylabel('QA')
xlabel('ncount');
title(['r=-0.1, p=0.08']);
idx = find(ncount>0);
[r,  p] = corr(qa(idx), log(ncount(idx)));