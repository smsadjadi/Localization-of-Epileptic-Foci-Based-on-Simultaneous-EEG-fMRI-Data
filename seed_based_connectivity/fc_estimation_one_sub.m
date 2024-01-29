close all; clear; clc;

%% Seed-Based Functional Connectivity

SubDataPath = 'D:\MEGA\Data\Sub_01\';
fMRIpath = [SubDataPath,'Outputs\M123outputs.feat\gunzip\'];

BOLD.hdr = spm_vol([fMRIpath,'filtered_func_data.nii']);
BOLD.savepath = [SubDataPath,'Outputs\M3'];
BOLD.img = spm_read_vols(BOLD.hdr);
BOLD.ID = 1;

mask = spm_read_vols(spm_vol([fMRIpath,'mask.nii']));

BOLDresponse.hdr = spm_vol([fMRIpath,'rendered_thresh_zstat1.nii']);
BOLDresponse.img = spm_read_vols(BOLDresponse.hdr);
[c1,c2,c3] = ind2sub(size(BOLDresponse.img),find(BOLDresponse.img>4));
c4 = BOLDresponse.img(BOLDresponse.img>4); BOLDcube = [c1,c2,c3,c4];
c5 = maxk(c4,10); c6 = ismember(BOLDcube, c5);
[c7,~] = find(c6); BOLDresponse.maxBOLDcube = BOLDcube(c7,:);
BOLDresponse.maxBOLDmni = cube2mni(BOLDresponse.maxBOLDcube(:,1:3),BOLDresponse.hdr(1));
BOLDresponse.maxBOLDmni(:,4) = BOLDresponse.maxBOLDcube(:,4);
clear('c1','c2','c3','c4','c5','c6','c7','BOLDcube')

for s = 1:10
    seed(s).name = ['maxBOLD',num2str(s)];
    seed(s).coordinates = BOLDresponse.maxBOLDmni(s,1:3);
    seed(s).connum = 27;
end

voxel_wise_connectivity(BOLD,seed,mask);

% load('D:\MEGA\Data\Sub_01\Outputs\M3\seed1_timecourse.txt')
% plot(seed1_timecourse); grid on

%% Cluster Mask Extraction
% 
% Zthresh = 1.7;
% Zmap = zeros(78,78,40);
% for s = 1:10
%     zmap(s).seed = spm_read_vols(spm_vol([BOLD.savepath,filesep,'seed',num2str(s),'FCzmap.nii']));
%     zmap(s).seed(zmap(s).seed<Zthresh) = 0;
%     Zmap = Zmap + zmap(s).seed;
% end
% Zmap(Zmap~=0)= 1;
% 
% Z2map = zeros(size(Zmap));
% for i = 2:76
%     for j = 2:76
%         for k = 2:38
%             cube = Zmap(i-1:i+1,j-1:j+1,k-1:k+1);
%             cube(1,1,1) = 0; cube(end,1,1) = 0; cube(1,end,1) = 0; cube(1,1,end) = 0;
%             cube(1,end,end) = 0; cube(end,1,end) = 0; cube(end,end,1) = 0; cube(end,end,end) = 0;
%             if sum(cube(:))>0
%                 Z2map(i-1:i+1,j-1:j+1,k-1:k+1)=1;
%             end
%         end
%     end
% end
% 
% Z3map = zeros(size(Z2map));
% for i = 3:75
%     for j = 3:75
%         for k = 3:37
%             cube = Z2map(i-2:i+2,j-2:j+2,k-2:k+2);
% %             cube(1,1,1) = 0; cube(end,1,1) = 0; cube(1,end,1) = 0; cube(1,1,end) = 0;
% %             cube(1,end,end) = 0; cube(end,1,end) = 0; cube(end,end,1) = 0; cube(end,end,end) = 0;
%             if sum(cube(:))==sum(sum(sum(ones(size(cube)))))
%                 Z3map(i,j,k)=1;
%             end
%         end
%     end
% end
% 
% hdr = BOLD.hdr(1);
% hdr.fname = [BOLD.savepath,filesep,'FCzmap.nii'];
% hdr.dt = [16,0];
% hdr = spm_write_vol(hdr,Z3map);