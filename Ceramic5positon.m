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
[s1_ss, s2_ss, s3_ss] = size(mtype_cube_subsamp);
model = mtype_cube_subsamp; air_id = -1;
sliceomatic(model)
% file_name = 'scat_fib_1.csv';
% create_csv(file_name,model,s1_ss,s2_ss,s3_ss,air_id,mtype_cube_subsamp,pval_cube_subsamp);

figure; colormap(gray);
contourf(model(:,:,floor(s3_ss/2))); 
% model_pval = pval_cube_subsamp;
figure; contourf((0:s2_ss-1)*.1,(0:s1_ss-1)*.1,model(:,:,floor(s3_ss/2)),'LineStyle','none');
% colormap(flipud(gray));brighten(.4);
xlabel('Distance (cm)','FontSize',14); ylabel('Distance (cm)','FontSize',14);
hcb = colorbar; 
set(hcb,'YTick',[0,1,2,3,4,5,6,7,8],'YTickLabel',{'Air','Skin','Gland-1','Gland-2','Gland-3','Fat-1','Fat-2','Fat-3','Muscle'},'FontSize',14);

% create temperature model
tumor_on = 0; tumor_depth = 8;  tumor_radius = 10; Tambient = 27; Tart = 37; 
% muscle_wall = 154; skin_start = 0; % can't remember why I put this in
[T_3d_nom,tissue_3d_nom] = gen_breast_therm_model(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start);
figure; contourf(T_3d_nom(:,:,floor(s3_ss/2)));
colorbar;
xlabel('Distance (mm)','FontSize',14); ylabel('Distance (mm)','FontSize',14);
title('Normal Temperature Profile (\circC)','FontSize',14);

% Generate temperature anomalies with radius = 10
tumor_on = 1; tum_y_cen = 90; tum_z_cen = floor(s3_ss/2);

tumor_radius = 10; tumor_depth = 10; tum_x_cen = 20 + tumor_depth;% tumor dept is 1cm
[T_3d_abn1,tissue_3d_abn1] = gen_breast_therm_model(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start,tum_x_cen,tum_y_cen,tum_z_cen);
tumor_radius = 10; tumor_depth = 20; tum_x_cen = 20 + tumor_depth;% tumor depth is 2cm
[T_3d_abn2,tissue_3d_abn2] = gen_breast_therm_model(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start,tum_x_cen,tum_y_cen,tum_z_cen);
tumor_radius = 10; tumor_depth = 30; tum_x_cen = 20 + tumor_depth;% tumor depth is 3cm
[T_3d_abn3,tissue_3d_abn3] = gen_breast_therm_model(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start,tum_x_cen,tum_y_cen,tum_z_cen);
tumor_radius = 10; tumor_depth = 40; tum_x_cen = 20 + tumor_depth;% tumor depth is 4cm
[T_3d_abn4,tissue_3d_abn4] = gen_breast_therm_model(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start,tum_x_cen,tum_y_cen,tum_z_cen);


figure; contourf(T_3d_abn1(:,:,floor(s3_ss/2)));
xlabel('Distance (mm)','FontSize',14); ylabel('Distance (mm)','FontSize',14);
hcb = colorbar; 
set(hcb,'YTick',[0,1,2,3,4,5,6,7,8],'YTickLabel',{'Air','Skin','Gland-1','Gland-2','Gland-3','Fat-1','Fat-2','Fat-3','Muscle'},'FontSize',14);
title('Temperature Profile of 1cm Deep Tumor (\circC)','FontSize',14);

figure; contourf(T_3d_abn2(:,:,floor(s3_ss/2)));
xlabel('Distance (mm)','FontSize',14); ylabel('Distance (mm)','FontSize',14);
hcb = colorbar; 
set(hcb,'YTick',[0,1,2,3,4,5,6,7,8],'YTickLabel',{'Air','Skin','Gland-1','Gland-2','Gland-3','Fat-1','Fat-2','Fat-3','Muscle'},'FontSize',14);
title('Temperature Profile of 2cm Deep Tumor (\circC)','FontSize',14);

figure; contourf(T_3d_abn3(:,:,floor(s3_ss/2)));
xlabel('Distance (mm)','FontSize',14); ylabel('Distance (mm)','FontSize',14);
hcb = colorbar; 
set(hcb,'YTick',[0,1,2,3,4,5,6,7,8],'YTickLabel',{'Air','Skin','Gland-1','Gland-2','Gland-3','Fat-1','Fat-2','Fat-3','Muscle'},'FontSize',14);
title('Temperature Profile of 3cm Deep Tumor (\circC)','FontSize',14);

figure; 
contourf(T_3d_abn4(:,:,floor(s3_ss/2)));
xlabel('Distance (mm)','FontSize',14); ylabel('Distance (mm)','FontSize',14);
hcb = colorbar; 
set(hcb,'YTick',[0,1,2,3,4,5,6,7,8],'YTickLabel',{'Air','Skin','Gland-1','Gland-2','Gland-3','Fat-1','Fat-2','Fat-3','Muscle'},'FontSize',14);
title('Temperature Profile of 4cm Deep Tumor (\circC)','FontSize',14);


% % % print models to csv file for xfdtd simulations
% % f=1;
% % create_csvfile(tissue_3d_nom,s1_ss,s2_ss,s3_ss,air_id,'scat_fib_nom.csv',pval_cube_subsamp,mtype_cube_subsamp,f);
% % create_csvfile(tissue_3d_abn1,s1_ss,s2_ss,s3_ss,air_id,'scat_fib_abn1.csv',pval_cube_subsamp,mtype_cube_subsamp,f);
% % create_csvfile(tissue_3d_abn2,s1_ss,s2_ss,s3_ss,air_id,'scat_fib_abn2.csv',pval_cube_subsamp,mtype_cube_subsamp,f);
% % create_csvfile(tissue_3d_abn3,s1_ss,s2_ss,s3_ss,air_id,'scat_fib_abn3.csv',pval_cube_subsamp,mtype_cube_subsamp,f);
% % create_csvfile(tissue_3d_abn4,s1_ss,s2_ss,s3_ss,air_id,'scat_fib_abn4.csv',pval_cube_subsamp,mtype_cube_subsamp,f)

% plot 1d temperatures as is
h_1d_temp = figure; 
hold on;
plot(T_3d_nom(:,90,80),'k','LineWidth',1.5);
plot(T_3d_abn1(:,90,80),'r','LineWidth',1.5);
plot(T_3d_abn2(:,90,80),'y','LineWidth',1.5);
plot(T_3d_abn3(:,90,80),'b','LineWidth',1.5);
plot(T_3d_abn4(:,90,80),'m','LineWidth',1.5);
xlabel('Distance (mm)','FontSize',14); ylabel('Temperature ( \circ C)','FontSize',14);
habn_temp_diff_legend = legend('normal','Tumor 1cm deep','Tumor 2cm deep','Tumor 3cm deep','Tumor 4cm','Location','Southeast'); 
set(habn_temp_diff_legend,'FontSize',8);
hold off;

% Plot 1d temperature difference between nominal and abnormal
h_1d_temp_diff = figure; hold on;
plot(T_3d_abn1(20:100,90,80)-T_3d_nom(20:100,90,80),'k','LineWidth',1);
plot(T_3d_abn2(20:100,90,80)-T_3d_nom(20:100,90,80),'k--','LineWidth',1);
plot(T_3d_abn3(20:100,90,80)-T_3d_nom(20:100,90,80),'k:','LineWidth',1);
plot(T_3d_abn4(20:100,90,80)-T_3d_nom(20:100,90,80),'k-*','LineWidth',1);
xlabel('Distance (mm)','FontSize',14); ylabel('Temperature Difference (\circC)','FontSize',14);
habn_temp_diff_legend = legend('Tumor 1cm deep','Tumor 2cm deep','Tumor 3cm deep','Tumor 4cm deep','Location','NorthEast'); 
set(habn_temp_diff_legend,'FontSize',8);
axis([0 100 0 3.8]);
hold off;
title('Temperature Increase due to Tumors at Different Depths','FontSize',14);

model = model(18:end,:,:);
[s1_ss, s2_ss, s3_ss] = size(model);
load scat_fib_nom_22ghz_ceramic_center1.mat;
load scat_fib_nom_22ghz_ceramic_right1.mat;
load scat_fib_nom_22ghz_ceramic_left1.mat;
load scat_fib_nom_22ghz_ceramic_top1.mat;
load scat_fib_nom_22ghz_ceramic_bottom1.mat;


% Read WFs into a cube
wf_cube_2_2ghz_ceramic_center = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_22ghz_ceramic_center1);
wf_cube_2_2ghz_ceramic_right = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_22ghz_ceramic_right1);
wf_cube_2_2ghz_ceramic_left = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_22ghz_ceramic_left1);
wf_cube_2_2ghz_ceramic_top = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_22ghz_ceramic_top1);
wf_cube_2_2ghz_ceramic_bottom = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_22ghz_ceramic_bottom1);

% eliminate extra data
wf_cube_2_2ghz_ceramic_center = wf_cube_2_2ghz_ceramic_center(:,1:s2_ss,:);
wf_cube_2_2ghz_ceramic_right = wf_cube_2_2ghz_ceramic_right(:,1:s2_ss,:);
wf_cube_2_2ghz_ceramic_left = wf_cube_2_2ghz_ceramic_left(:,1:s2_ss,:);
wf_cube_2_2ghz_ceramic_top = wf_cube_2_2ghz_ceramic_top(:,1:s2_ss,:);
wf_cube_2_2ghz_ceramic_bottom = wf_cube_2_2ghz_ceramic_bottom(:,1:s2_ss,:);

% find norms of wfs and plot cube results
% convert cube to 1d
wf_vec_2_2ghz_ceramic_center = convert_3d_to_1d(wf_cube_2_2ghz_ceramic_center,s1_ss,s2_ss,s3_ss);
wf_vec_2_2ghz_ceramic_right = convert_3d_to_1d(wf_cube_2_2ghz_ceramic_right,s1_ss,s2_ss,s3_ss);
wf_vec_2_2ghz_ceramic_left = convert_3d_to_1d(wf_cube_2_2ghz_ceramic_left,s1_ss,s2_ss,s3_ss);
wf_vec_2_2ghz_ceramic_top = convert_3d_to_1d(wf_cube_2_2ghz_ceramic_top,s1_ss,s2_ss,s3_ss);
wf_vec_2_2ghz_ceramic_bottom = convert_3d_to_1d(wf_cube_2_2ghz_ceramic_bottom,s1_ss,s2_ss,s3_ss);

% sum and find the norm of each wf
wf_vec_2_2ghz_norm_ceramic_center = wf_vec_2_2ghz_ceramic_center/sum(wf_vec_2_2ghz_ceramic_center);
wf_vec_2_2ghz_norm_ceramic_right = wf_vec_2_2ghz_ceramic_center/sum(wf_vec_2_2ghz_ceramic_right);
wf_vec_2_2ghz_norm_ceramic_left = wf_vec_2_2ghz_ceramic_center/sum(wf_vec_2_2ghz_ceramic_left);
wf_vec_2_2ghz_norm_ceramic_top = wf_vec_2_2ghz_ceramic_center/sum(wf_vec_2_2ghz_ceramic_top);
wf_vec_2_2ghz_norm_ceramic_bottom = wf_vec_2_2ghz_ceramic_center/sum(wf_vec_2_2ghz_ceramic_bottom);

wf_cube_2_2ghz_norm_ceramic_center = convert_1d_to_3d(wf_vec_2_2ghz_norm_ceramic_center,s1_ss,s2_ss,s3_ss);% 1.5ghz
wf_cube_2_2ghz_norm_ceramic_right = convert_1d_to_3d(wf_vec_2_2ghz_norm_ceramic_right,s1_ss,s2_ss,s3_ss);% 1.5ghz
wf_cube_2_2ghz_norm_ceramic_left = convert_1d_to_3d(wf_vec_2_2ghz_norm_ceramic_left,s1_ss,s2_ss,s3_ss);% 1.5ghz
wf_cube_2_2ghz_norm_ceramic_top = convert_1d_to_3d(wf_vec_2_2ghz_norm_ceramic_top,s1_ss,s2_ss,s3_ss);% 1.5ghz
wf_cube_2_2ghz_norm_ceramic_bottom = convert_1d_to_3d(wf_vec_2_2ghz_norm_ceramic_bottom,s1_ss,s2_ss,s3_ss);% 1.5ghz

% Plot Weighting Functions With Normalization
figure
semilogy(wf_cube_2_2ghz_norm_ceramic_center(:,floor(s2_ss/3),floor(s3_ss/3)),'k','LineWidth',1.5);hold on;
semilogy(wf_cube_2_2ghz_norm_ceramic_right(:,floor(s2_ss/3),floor(s3_ss/3)),'g','LineWidth',1.5);
semilogy(wf_cube_2_2ghz_norm_ceramic_left(:,floor(s2_ss/3),floor(s3_ss/3)),'b','LineWidth',1.5); 
semilogy(wf_cube_2_2ghz_norm_ceramic_top(:,floor(s2_ss/3),floor(s3_ss/3)),'m','LineWidth',1.5); 
semilogy(wf_cube_2_2ghz_norm_ceramic_bottom(:,floor(s2_ss/3),floor(s3_ss/3)),'r','LineWidth',1.5);

xlabel('Depth (mm)');ylabel('Weighting Functions');
title('Weight functions with normalization');
habn_wf_legend = legend('center','right','left','top','bottom','Location','Northeast'); 
set(habn_wf_legend,'FontSize',8);
hold off;

WF = [wf_vec_2_2ghz_norm_ceramic_center,wf_vec_2_2ghz_norm_ceramic_right,wf_vec_2_2ghz_norm_ceramic_left,wf_vec_2_2ghz_norm_ceramic_top,wf_vec_2_2ghz_norm_ceramic_bottom];
% Calculate Nominal Brightness Temperature, TB_nominal
Tvec_nominal = convert_3d_to_1d(T_3d_nom,s1_ss,s2_ss,s3_ss);
TB_nominal = (WF'*Tvec_nominal)';
% Calculate Brightness Temperature with Tumors
Tvec_abn1 = convert_3d_to_1d(T_3d_abn1,s1_ss,s2_ss,s3_ss);
Tvec_abn2 = convert_3d_to_1d(T_3d_abn2,s1_ss,s2_ss,s3_ss);
Tvec_abn3 = convert_3d_to_1d(T_3d_abn3,s1_ss,s2_ss,s3_ss);
Tvec_abn4 = convert_3d_to_1d(T_3d_abn4,s1_ss,s2_ss,s3_ss);

TB_abn1 = (WF'*Tvec_abn1)';
TB_abn2 = (WF'*Tvec_abn2)';
TB_abn3 = (WF'*Tvec_abn3)';
TB_abn4 = (WF'*Tvec_abn4)';

% Plot Brightness Temperature difference between Normal and Abnormal Breast
TB_diff1 = TB_abn1-TB_nominal;
TB_diff2 = TB_abn2-TB_nominal;
TB_diff3 = TB_abn3-TB_nominal;
TB_diff4 = TB_abn4-TB_nominal;
Y = [TB_diff1;TB_diff2;TB_diff3;TB_diff4];
TB_diff_plot = figure;
%plot(TB_diff1,'k','LineWidth',1); hold on;
%plot(TB_diff2,'k:','LineWidth',1); hold on;
%plot(TB_diff3,'k--','LineWidth',1); hold on;
%plot(TB_diff4,'r','LineWidth',1); hold on;
%xlabel('Different Frequencies and Positions','FontSize',14); ylabel('Temperature Difference (\circC)','FontSize',14);
%habn_temp_diff_legend = legend('Tumor 1cm deep','Tumor 2cm deep','Tumor 3cm deep','Tumor 4cm deep','Location','Northeast'); 
%set(habn_temp_diff_legend,'FontSize',8);
%title('Brightness Temperature Difference','FontSize',14);
%axis([1 4 0 0.1]);
%hold off;
plot(Y);
xlabel('Tumor Depth','FontSize',14); ylabel('Temperature Difference (\circC)','FontSize',14);
habn_temp_diff_legend = legend('center','right','left','','top','bottom','Northeast'); 
set(habn_temp_diff_legend,'FontSize',8);
title('Brightness Temperature Difference','FontSize',14);
%axis([1 4 0 0.1]);
%hold off;
% calculate SNR

TB_diff1_avg = sum(TB_diff1)/5;
TB_diff2_avg = sum(TB_diff2)/5;
TB_diff3_avg = sum(TB_diff3)/5;
TB_diff4_avg = sum(TB_diff4)/5;

Noise = 0.1;
SNR1 =  20*log(((TB_diff1_avg^2)/5)/((Noise/5)^2));
SNR2 =  20*log(((TB_diff2_avg^2)/5)/((Noise/5)^2));
SNR3 =  20*log(((TB_diff3_avg^2)/5)/((Noise/5)^2));
SNR4 =  20*log(((TB_diff4_avg^2)/5)/((Noise/5)^2));

