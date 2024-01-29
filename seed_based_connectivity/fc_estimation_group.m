clear;
clc;
firstlevel_path = 'D:\VMWare_shared\MNI_research\selected_first_result';
functional_path = 'D:\VMWare_shared\MNI_research\selected_functional';
[num,txt,raw] = xlsread('D:\VMWare_shared\MNI_research\code\FC_IED\selection.xls');
subject = txt;
% Type of study (1,2,3 in the article)
study = num(:,1);
% This is actually not used, can be ignored
feature = num(:,2);
% MNI coordinate of the maximal BOLD response 
mnicoor = num(:,3:5);
transformcoor = zeros(length(study),3);

for s = 1:length(subject)
    prename = subject{s,1}(1:5);
    subject_folder = dir([functional_path,filesep,prename,'*']);
    sub_functional_path = fullfile(functional_path,subject_folder(1).name,'functional_nii');
    firstlevel_folder = dir([firstlevel_path,filesep,prename,'*']);
    sub_firstlevel_path = fullfile(firstlevel_path,firstlevel_folder(1).name);
    cd(sub_functional_path);
    
    fc_folder = fullfile(functional_path,subject_folder(1).name,['sm_fcMRI_study_',sprintf('%02d',study(s))]);
    if ~exist(fc_folder)
        mkdir(fc_folder);
    end
    % seed.name: region name
    % seed.coordinates: (x,y,z), MNI space
    % seed.connum: number of neighbourhood, 27 or 81 voxels
    seed.name = ['study_',sprintf('%02d',study(s))];
    % This is the file of the normalized T1 image (to MNI space) of each subject
    native_r = fullfile(sub_firstlevel_path,'T1_native_r.nii');
    if exist(native_r)
        hdr_r = spm_vol(native_r);
        % This is the file of the original T1 image (in subject space) of each subject
        hdr_native = spm_vol(fullfile(sub_firstlevel_path,'T1_native.nii'));
        tranform_mat = inv(hdr_native.mat)*hdr_r.mat;
        seed.coordinates = (tranform_mat*[mnicoor(s,:),1]')';
        seed.coordinates(4) = [];
    else
        seed.coordinates = mnicoor(s,:);
    end
    seed.connum = 27;
    transformcoor(s,:) = seed.coordinates;
    
    % This is all the files of the processed (slice timing, head motion correction, temporal and spatial smoothing, temporal filtering) fmri images of each subject
    functional_file = dir([sub_functional_path,filesep,'smoothed_filtered',filesep,'*mri.nii']);
    for f = 1:length(functional_file)
        fprintf('calculating study %d scan %d \n',s,f);
        % subject.hdr: the header file of the BOLD images
        % subject.ID;
        % subject.img: the BOLD images of subjects
        % subject.savepath:
        BOLD.hdr = spm_vol([sub_functional_path,filesep,'smoothed_filtered',filesep,functional_file(f).name]);
        BOLD.ID = [prename,'_scan_',sprintf('%02d',f)];
        BOLD.img = spm_read_vols(BOLD.hdr);
        BOLD.savepath = [fc_folder,filesep];
        VoxelWiseConnectivityCalculate(BOLD,seed);
    end
end