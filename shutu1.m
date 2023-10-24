%%
clear;
clc;

I_target = imread('D:\�ļ�\����ͼ����\ͼ��ʶ��������Ŀ\target2.png');          %����ͼƬ
I_template = imread('D:\�ļ�\����ͼ����\ͼ��ʶ��������Ŀ\template2.png');

figure
subplot(121)
imshow(I_target);
title('ʵ��ͼƬ')
subplot(122)
imshow(I_template);
title('��Ƭ')

v_target = rgb2gray(I_target);           % ͼ��ת��      
v_template = rgb2gray(I_template);


j_target = medfilt2(v_target,[3,3]);      %3��3��ֵ�˲���ȥ��
j_template = medfilt2(v_template,[3,3]); 

% figure
% subplot(121)
% imshow(j_target);
% title('��ֵ�˲�ȥ����ʵ��ͼƬ')
% subplot(122)
% imshow(j_template);
% title('��ֵ�˲�ȥ������Ƭ')

j_target = imadjust(j_target, [], [], 0.75);
j_template = imadjust(j_template, [], [], 0.75);

% figure
% subplot(121)
% imshow(j_target);
% title('٤��任���ʵ��ͼƬ')
% subplot(122)
% imshow(j_template);
% title('٤��任�����Ƭ')

[m,n]=size(j_target);             %���ߴ�
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
T = 0.978;   %������ص���ֵ
[iMaxPos,jMaxPos]=find( result>=T); 
N = length(iMaxPos);

%%
I_target2 = I_target;
wSize = [n0,m0];
lineSize = 1;
color = [0,0,0];
for i = 1:N    
%     rectangle('position',[jMaxPos(i),iMaxPos(i),n0,n0],'edgecolor','w');      
%     plot(jMaxPos(i),iMaxPos(i),'*');              %������ص�  
%     plot([jMaxPos(i),jMaxPos(i)+n0-1],[iMaxPos(i),iMaxPos(i)]);     %�þ��ο��ǳ�ƥ������
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
disp(['��������������������������������']);
disp(['оƬʶ��ɹ���']);
disp(['ʶ��ɹ���оƬ����Ϊ' ,num2str(num)]);
disp(['��������������������������������']);

%%
Ilabel = bwlabel(Ibw2,8);
stat = regionprops(Ilabel,'centroid');  %��������
stat(1) = [];

m = length(stat);
x1 = zeros(m,1);
x2 = zeros(m,1);
y1 = zeros(m,1);
y2 = zeros(m,1);

figure
imshow(I_target); 
title('ʶ��ɹ���оƬ')
hold on;

for x = 1: numel(stat)
    plot(stat(x).Centroid(1),stat(x).Centroid(2),'g+');
end

%��������
for i = 1: numel(stat)
    x1(i) = round(stat(i).Centroid(1) - n0/2);
    x2(i) = round(stat(i).Centroid(1) + n0/2);
    y1(i) = round(stat(i).Centroid(2) - m0/2);
    y2(i) = round(stat(i).Centroid(2) + m0/2);
end

%���оƬ����
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
xlswrite('C:\Users\FENGYIQING\Desktop\ͼ��ʶ��������Ŀ\�ڶ������С������.xlsx',xuhao,'Sheet1','A2');
xlswrite('C:\Users\FENGYIQING\Desktop\ͼ��ʶ��������Ŀ\�ڶ������С������.xlsx',zuoshang,'Sheet1','B2');
xlswrite('C:\Users\FENGYIQING\Desktop\ͼ��ʶ��������Ŀ\�ڶ������С������.xlsx',youshang,'Sheet1','C2');
xlswrite('C:\Users\FENGYIQING\Desktop\ͼ��ʶ��������Ŀ\�ڶ������С������.xlsx',zuoxia,'Sheet1','D2');
xlswrite('C:\Users\FENGYIQING\Desktop\ͼ��ʶ��������Ŀ\�ڶ������С������.xlsx',youxia,'Sheet1','E2');
xlswrite('C:\Users\FENGYIQING\Desktop\ͼ��ʶ��������Ŀ\�ڶ������С������.xlsx',zhongdian,'Sheet1','F2');

