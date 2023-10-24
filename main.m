%%
clear;
clc;

I_target = imread('target2.png');          %读入图片
I_template = imread('template2.png');

figure
subplot(121)
imshow(I_target);
title('实采图片')
subplot(122)
imshow(I_template);
title('样片')

v_target = rgb2gray(I_target);           % 图像转换      
v_template = rgb2gray(I_template);


j_target = medfilt2(v_target,[3,3]);      %3×3中值滤波器去噪
j_template = medfilt2(v_template,[3,3]); 

% figure
% subplot(121)
% imshow(j_target);
% title('中值滤波去噪后的实采图片')
% subplot(122)
% imshow(j_template);
% title('中值滤波去噪后的样片')

j_target = imadjust(j_target, [], [], 0.75);
j_template = imadjust(j_template, [], [], 0.75);

% figure
% subplot(121)
% imshow(j_target);
% title('伽马变换后的实采图片')
% subplot(122)
% imshow(j_template);
% title('伽马变换后的样片')

[m,n]=size(j_target);             %读尺寸
[m0,n0]=size(j_template);

result=zeros(m-m0+1,n-n0+1);
vec_sub = double( j_template(:) );
norm_sub = norm( vec_sub );

for i=1:1:m-m0+1
    for j=1:1:n-n0+1
        subMatr=j_target(i:i+m0-1,j:j+n0-1);
        vec=double( subMatr(:) );
        result(i,j)=vec'*vec_sub / (norm(vec)*norm_sub+eps);
    end
end

%%
T = 0.978;   %设置相关的阈值
[iMaxPos,jMaxPos]=find( result>=T); 
N = length(iMaxPos);

%%
I_target2 = I_target;
wSize = [n0,m0];
lineSize = 1;
color = [0,0,0];
for i = 1:N    
%     rectangle('position',[jMaxPos(i),iMaxPos(i),n0,n0],'edgecolor','w');      
%     plot(jMaxPos(i),iMaxPos(i),'*');              %绘制相关点  
%     plot([jMaxPos(i),jMaxPos(i)+n0-1],[iMaxPos(i),iMaxPos(i)]);     %用矩形框标记出匹配区域
%     plot([jMaxPos(i)+n0-1,jMaxPos(i)+n0-1],[iMaxPos(i),iMaxPos(i)+m0-1]);
%     plot([jMaxPos(i),jMaxPos(i)+n0-1],[iMaxPos(i)+m0-1,iMaxPos(i)+m0-1]);
%     plot([jMaxPos(i),jMaxPos(i)],[iMaxPos(i),iMaxPos(i)+m0-1]);
      pt = [jMaxPos(i),iMaxPos(i)];
      I_target2 = drawRect(I_target2, pt, wSize, lineSize, color);
end
figure
imshow(I_target2);
Ibw = imbinarize(I_target2,0.01);    

Ibw = im2uint8(Ibw);
figure
imshow(Ibw)

Ibw1 = rgb2hsv(Ibw);
Ibw2 = mat2gray(Ibw1(:,:,3));
[L,num] = bwlabel(Ibw2,8);
num = num - 1;
disp(['————————————————']);
disp(['芯片识别成功！']);
disp(['识别成功的芯片数量为' ,num2str(num)]);
disp(['————————————————']);

%%
Ilabel = bwlabel(Ibw2,8);
stat = regionprops(Ilabel,'centroid');  %形心坐标
stat(1) = [];

m = length(stat);
x1 = zeros(m,1);
x2 = zeros(m,1);
y1 = zeros(m,1);
y2 = zeros(m,1);

figure
imshow(I_target); 
title('识别成功的芯片')
hold on;

for x = 1: numel(stat)
    plot(stat(x).Centroid(1),stat(x).Centroid(2),'g+');
end

%计算坐标
for i = 1: numel(stat)
    x1(i) = round(stat(i).Centroid(1) - n0/2);
    x2(i) = round(stat(i).Centroid(1) + n0/2);
    y1(i) = round(stat(i).Centroid(2) - m0/2);
    y2(i) = round(stat(i).Centroid(2) + m0/2);
end

%框出芯片区域
for i = 1: numel(stat)
    rectangle('position',[x1(i),y1(i),n0,m0],'edgecolor','r','LineWidth',1);  
end

%%
x1 = string(x1);
x2 = string(x2);
y1 = string(y1);
y2 = string(y2);

zuoshang = zeros(m,1);
zuoxia = zeros(m,1);
youshang = zeros(m,1);
youxia = zeros(m,1);
zhongdian = zeros(m,1);
zuoshang = string(zuoshang);
zuoxia = string(zuoxia);
youshang = string(youshang);
youxia = string(youxia);
zhongdian = string(zhongdian);

for i = 1: numel(stat)
    zuoshang(i) = strcat(x1(i),' , ',y1(i));
    zuoxia(i) = strcat(x1(i),' , ',y2(i));
    youshang(i) = strcat(x2(i),' , ',y1(i));
    youxia(i) = strcat(x2(i),' , ',y2(i));
    zhongdian(i) = strcat(string(stat(i).Centroid(1)),' , ',string(stat(i).Centroid(2)));
end
xuhao = 1:num;
xuhao = string(xuhao');

