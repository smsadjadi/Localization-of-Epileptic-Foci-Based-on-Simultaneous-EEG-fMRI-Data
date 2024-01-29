function second_level_test(result_path,zmap_file,ncontra,contrastvec,contrast_name,significance,extent)

% spm_get_defaults
global defaults
spm('defaults','fmri');
spm_jobman('initcfg');
% tail = '2nd_level';
% ncontra = 1;
% contrast_name = {'pre-NC','NC-pre','NC-post','post-NC'};
% contrastvec = {1,-1};

%%
jobs{1}.stats{1}.factorial_design.dir = cellstr(result_path);
jobs{1}.stats{1}.factorial_design.des.t1.scans = zmap_file;

jobs{1}.stats{1}.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
jobs{1}.stats{1}.factorial_design.masking.tm.tm_none = 1;
jobs{1}.stats{1}.factorial_design.masking.im = 0;
jobs{1}.stats{1}.factorial_design.masking.em = {''};
jobs{1}.stats{1}.factorial_design.globalc.g_omit = 1;
jobs{1}.stats{1}.factorial_design.globalm.gmsca.gmsca_no = 1;
jobs{1}.stats{1}.factorial_design.globalm.glonorm = 1;

%% 
jobs{1}.stats{2}.fmri_est.spmmat = cellstr(fullfile([result_path,'\SPM.mat']));
jobs{1}.stats{2}.fmri_est.method.Classical = 1;
% inputs = cell(0, 1);
% spm_jobman('serial', jobs, '', inputs{:});
spm_jobman('run',jobs);
clear jobs;
%% contrast manager
matlabbatch{1}.spm.stats.con.spmmat = cellstr(fullfile([result_path,'\SPM.mat']));
for con = 1:ncontra
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name = contrast_name;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.convec = contrastvec;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'none';
end
matlabbatch{1}.spm.stats.con.delete = 0;
%% results report
matlabbatch{2}.spm.stats.results.spmmat = cellstr(fullfile([result_path,'\SPM.mat']));
for con = 1:ncontra
    matlabbatch{2}.spm.stats.results.conspec(con).titlestr = contrast_name;
    matlabbatch{2}.spm.stats.results.conspec(con).contrasts = con;
    matlabbatch{2}.spm.stats.results.conspec(con).threshdesc = 'none';
    matlabbatch{2}.spm.stats.results.conspec(con).thresh = significance;
    matlabbatch{2}.spm.stats.results.conspec(con).extent = extent;
    matlabbatch{2}.spm.stats.results.conspec(con).conjunction = 1;
%     matlabbatch{2}.spm.stats.results.conspec(con).mask.image.name = cellstr(maskfile);
%     matlabbatch{2}.spm.stats.results.conspec(con).mask.image.mtype = 0;
end
matlabbatch{2}.spm.stats.results.units = 1;
matlabbatch{2}.spm.stats.results.export{1}.fig = true;
matlabbatch{2}.spm.stats.results.export{2}.binary.basename = contrast_name;
spm_jobman('run',matlabbatch);
% tail = 'spm_filtered';
% F = contrast_name{1};
% Vo = spm_write_filtered(xSPM.Z,xSPM.XYZ,xSPM.DIM,xSPM.M,tail,F);
clear matlabbatch;
end     