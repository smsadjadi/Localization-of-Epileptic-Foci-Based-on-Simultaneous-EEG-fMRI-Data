function [] = voxel_wise_connectivity(BOLD,seed,brainmask)

% Usage:
% fcmatrx = VoxelWiseConnectivityCalculate(4Ddata,seed,brainmask)
% BOLD:
% BOLD.hdr: the header file of the BOLD images
% BOLD.ID: ...
% BOLD.img: the BOLD images of subjects
% BOLD.savepath: ...
% seed:
% seed.name: region name
% seed.coordinates: (x,y,z), MNI space
% seed.connum: number of neighbourhood, 27 or 81 voxels
% brainmask: The binary mask of the brain volume corresponding to the BOLD
% volume size.

hdr = BOLD.hdr(1);
data = BOLD.img;
savepath = BOLD.savepath;
subjectID = BOLD.ID;
[x,y,z,L] = size(data);
seednum = length(seed);
for s = 1:seednum
    seed(s).cubecoor = mni2cube(seed(s).coordinates,hdr);
end
for f = 1:L
    temp = data(:,:,:,f);
    for s = 1:seednum
        if seed(s).connum == 27
            littlecube = temp(seed(s).cubecoor(1)-1:seed(s).cubecoor(1)+1,...
                              seed(s).cubecoor(2)-1:seed(s).cubecoor(2)+1,....
                              seed(s).cubecoor(3)-1:seed(s).cubecoor(3)+1);
        elseif seed(s).connum == 81
            littlecube = temp(seed(s).cubecoor(1)-2:seed(s).cubecoor(1)+2,...
                              seed(s).cubecoor(2)-2:seed(s).cubecoor(2)+2,....
                              seed(s).cubecoor(3)-2:seed(s).cubecoor(3)+2);
        else
            error('Connection number should be 27 or 81\n');
        end
        seed_timecourse(f,s) = mean(littlecube(:));
    end
end
clear temp;
for s = 1:seednum
    timecourse = seed_timecourse(:,s);
    save([savepath,filesep,'seed',num2str(s),'_timecourse.txt'],'timecourse','-ASCII');
end
fcmatrix = zeros(x,y,z);
disp(['Calculating correlation ...']);
for i = 1:x   
    for j = 1:y
        for k = 1:z
            for s = 1:seednum
                temp(:,1) = data(i,j,k,:);
                if length(find(temp==0)) == length(temp)
                    seed(s).fcmatrix(i,j,k) = 0;
                    seed(s).zmap(i,j,k) = 0;
                else
                    seed(s).fcmatrix(i,j,k) = corr(seed_timecourse(:,s),temp);
                    seed(s).zmap(i,j,k) = 0.5*log((1+seed(s).fcmatrix(i,j,k))/(1-seed(s).fcmatrix(i,j,k)));
                end
            end
        end
    end
    disp(['slice ',num2str(i),'/',num2str(x),' is being calculated...']);     
end

disp('Saving Z-map...');
for s = 1:seednum
    seed(s).zmap = seed(s).zmap.*brainmask;
    % str = strcat(savepath,subjectID,'_',num2str(seed(s).coordinates(1)),'_',num2str(seed(s).coordinates(2)),'_',num2str(seed(s).coordinates(3)));
    hdr.fname = [savepath,filesep,'seed',num2str(s),'FCzmap.nii'];
    hdr.dt = [16,0];
    vector = isnan(seed(s).zmap);
    seed(s).zmap(vector==1) = 0;
    % seed(s).zmap(seed(s).zmap<1.7) = 0;
    hdr = spm_write_vol(hdr,seed(s).zmap);    
end
