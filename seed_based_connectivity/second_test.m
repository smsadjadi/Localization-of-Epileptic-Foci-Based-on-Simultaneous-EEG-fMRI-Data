clear;
clc;
datapath = 'D:\VMWare_shared\MNI_research\selected_functional';
subject = dir([datapath,filesep,'sub*']);

for s = 14%:length(subject)
    subfolder = fullfile(datapath,subject(s).name);
    study = dir([subfolder,filesep,'random_fcMRI_study*']);
    for f = [1 3]%:length(study)
        subsubfolder = fullfile(subfolder,study(f).name);
        zmap_str = spm_select('list', subsubfolder, '^*\zmap.nii');
        zmap_file = cell(size(zmap_str,1),1);
        for z = 1:size(zmap_str,1)
            zmap_file{z,1} = [subsubfolder,filesep,zmap_str(z,:),',1'];
        end

        contrastvec = {1;-1};
        ncontra = length(contrastvec);
        contrast_name = {'positive';'negative'};
        for re = 1:ncontra
            result_path = fullfile(subsubfolder,contrast_name{re,1});
            if ~exist(result_path)
                mkdir(result_path);
            end
            significance = 0.001;
            extent = 10;
            second_level_test(result_path,zmap_file,1,contrastvec{re,1},...
                contrast_name{re,1},significance,extent);
            F = 'result_mask';
            descrip = sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',xSPM.STAT,xSPM.u,xSPM.k);
            Vo = spm_write_filtered(xSPM.Z,xSPM.XYZ,xSPM.DIM,xSPM.M,descrip,F);
        end
    end
end