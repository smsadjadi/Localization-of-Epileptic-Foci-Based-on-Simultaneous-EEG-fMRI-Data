function [outpoints] = cube2mni(inpoints,V)
% Converts coordinates from MNI to cube brain
% FORMAT outpoints = mni2cube(inpoints)
% Where inpoints is N by 3 or 3 by N matrix of coordinates
%  (N being the number of points)
% outpoints is the coordinate matrix with cube points
% Ling-Li Zeng 9/2/11
% o is the origin of the brain, which is the cube coordinate corresponding
% to [0,0,0] in MNI coordinate.
% inpoints: input MNI coordinates, e.g. [3 -4 5]
% o: MNI origin in voxel coordinate, i.e. the voxel coordinates
% corresponding to where MNI = [0 0 0]
% V: image structure, reading using spm_vol function
% vs: voxel size, e.g. [3 3 3]
dimdim = find(size(inpoints) == 3);
if isempty(dimdim)
    error('input must be a N by 3 or 3 by N matrix');
elseif ~isempty(find(inpoints<0))
    error('no negative elements');
end
mat = V.mat;
[x,y] = size(inpoints);
outpoints = zeros(x,y+1);
for i = 1:1:size(inpoints,1)
    outpoints(i,:) = round(mat*[inpoints(i,:),1]')';
end
outpoints(:,4) = [];