clc
clear
addpath('F:/OneDrive/library/viggin-domain-adaptation-toolbox')
load('E:/projects/TPD/data/TPD_5percentage_mice_imputate.mat')
load('E:/projects/TPD/data/TPD_5percentage_mice_imputate_labels.mat')
maSrc=true(1,580);
maSrc(400:580)=false;

param = []; param.kerName = 'lin';param.bSstca = 0;
param.mu = 1;param.m = 218;param.gamma = .1;param.lambda = 0;
[newTPD,transMdl] = ftTrans_tca(TPD,maSrc,Labels(:,2),Labels(1,:),param);

dlmwrite('E:/projects/TPD/results/TPD_TCA_transformed.txt',newTPD,'delimiter','\t','precision',4);