clear;
load('heatingRes.mat');   %改mat名字
A = kdiff.sensor_data;
idx = 1;
% for i = 19:60:319      % 改4张图时间
figure;
for i = 5:5:80
    T1{idx} = reshape(A(:,i),490,490);
    subplot(4,4,idx);
    imagesc([0:22:220],[0:22:220],T1{idx});axis xy;
    idx = idx + 1;
end
figure;
subplot(2,3,1);
imagesc([0:22:220],[0:22:220],T1{1});axis xy;

subplot(2,3,2);
imagesc([0:22:220],[0:22:220],T1{2});axis xy;

subplot(2,3,3);
imagesc([0:22:220],[0:22:220],T1{3});axis xy;

subplot(2,3,4);
imagesc([0:22:220],[0:22:220],T1{4});axis xy;

subplot(2,3,5);
imagesc([0:22:220],[0:22:220],T1{5});axis xy;