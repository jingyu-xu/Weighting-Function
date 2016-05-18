fid = fopen('breastinfo_simple.txt','r');
breastID = str2double(fgetl(fid));
s1 = str2double(fgetl(fid));
s2 = str2double(fgetl(fid));
s3 = str2double(fgetl(fid));
class = str2double(fgetl(fid));
fclose(fid);

load mtype.mat;
load pval.mat;

muscle_wall = 153;
skin_start = 138;
z = [];
% Convert vector into cube
mtype_cube = zeros(s1,s2,s3); % each voxel is .5mmx.5mmx.5mm
pval_cube = zeros(s1,s2,s3);
cur_pos = 1;
for k=1:s3
    for j=1:s2
        for i= 1:s1
            mtype_cube(i,j,k) = mtype(cur_pos);
            pval_cube(i,j,k) = pval(cur_pos);
            cur_pos = cur_pos + 1;
        end 
    end
end

% subsample cubes in order to solve sparse matrix
s1_ss = floor(s1/2); % voxels are now 1mmx1mmx1mm
s2_ss = floor(s2/2);
s3_ss = floor(s3/2);
xi = 1; yi = 1; zi = 1;
mtype_cube_subsamp = zeros(s1_ss,s2_ss,s3_ss);
pval_cube_subsamp = zeros(s1_ss,s2_ss,s3_ss);
for z=1:2:s3-1
    for y = 1:2:s2
        for x = 1:2:s1
            mid = mtype_cube(x,y,z);
            pid = pval_cube(x,y,z);
            mtype_cube_subsamp(xi,yi,zi) = mid;
            pval_cube_subsamp(xi,yi,zi) = pid;
            xi = xi+1;
        end
        xi = 1;
        yi = yi + 1;
    end
    yi = 1;
    zi = zi + 1;
end
% some voxels not converted to muscle during subsampling
% so do that now
% still need to figure out how to get pval converted for fdtd
for x=1:s1_ss
    for y=1:s2_ss
        for z=1:s3_ss
            if x > 153
               mtype_cube_subsamp(x,y,z) = -4;
                pval_cube_subsamp(x,y,z) = 1;
            end
        end
    end
end

% save mtype_cube_subsamp; save pval_cube_subsamp.mat;
% load mtype_cube_subsamp.mat;
% [s1_ss, s2_ss, s3_ss] = size(mtype_cube_subsamp);
model = mtype_cube_subsamp; air_id = -1;
sliceomatic(model)
model = model(18:end,:,:);
[s1_ss, s2_ss, s3_ss] = size(model);

load scat_fib_nom_22ghz_ceramic_bottom1.mat;

% Read WFs into a cube
wf_cube_2_2ghz_ceramic_center = read_wf(s1_ss,s2_ss+1,s3_ss, scat_fib_nom_22ghz_ceramic_bottom1);
% eliminate extra data
wf_cube_2_2ghz_ceramic_center = wf_cube_2_2ghz_ceramic_center(:,1:s2_ss,:);
sliceomatic(wf_cube_2_2ghz_ceramic_center)
% find norms of wfs and plot cube results
% convert cube to 1d
wf_vec_2_2ghz_ceramic_center = convert_3d_to_1d(wf_cube_2_2ghz_ceramic_center,s1_ss,s2_ss,s3_ss);

% sum and find the norm of each wf
wf_vec_2_2ghz_norm_ceramic_center = wf_vec_2_2ghz_ceramic_center/sum(wf_vec_2_2ghz_ceramic_center);

wf_cube_2_2ghz_norm_ceramic_center = convert_1d_to_3d(wf_vec_2_2ghz_norm_ceramic_center,s1_ss,s2_ss,s3_ss);
for i = 1:s3_ss
    z(i) = wf_cube_2_2ghz_norm_ceramic_center(floor(s1_ss/3),floor(s2_ss/3),i);
end

% Plot Weighting Functions With Normalization
figure
semilogy((1:s3_ss), z,'k','LineWidth',1.5); hold on;
xlabel('Depth (mm)');ylabel('Weighting Functions');
title('3 GHz Weight functions from bottom');













