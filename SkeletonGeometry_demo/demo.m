clear all
clc
%% �����ͨ�������������ǵ�

P= load('C:\Users\Administrator\Desktop\���ܹ�������\��������ȡ\ffrReq.xyz');
info1 = load('C:\Users\Administrator\Desktop\���ܹ�������\��������ȡ\info.txt');
CoronaryInfo.PixelSpacing = info1(1:2);
CoronaryInfo.SpacingBetweenSlices = info1(3);

P(:,1) = round(P(:,1));
P(:,2) = round(P(:,2));
P(:,3) = round(P(:,3));
% �Ƿ���ĩ��
CoronaryInfo.CropEnds = 0;
minx =1;maxx = max(P(:,1))+5;
miny =1;maxy = max(P(:,2))+5;
minz =1;maxz = max(P(:,3))+5;
xvec = minx:1:maxx;
yvec = miny:1:maxy;
zvec = minz:1:maxz;
w = length(xvec);
l = length(yvec);
h = length(zvec);

ind = sub2ind([l,w,h],P(:,2),P(:,1),P(:,3));
%  ind = sub2ind([512,512,300],P(:,1),P(:,2),P(:,3));
CoronaryVol = false(l,w,h);
CoronaryVol(ind)=true;
 sum(sum(sum(CoronaryVol )))
 se = strel('sphere',1);
 CoronaryVol =imdilate(CoronaryVol, se);
% plot3(P(:,1),P(:,2),P(:,3),'g.','markersize',8,'linewidth',1)

CoronaryPara=[261 250 53 0 0;178 189 84 0 0];
%% CoronaryVol ������ά���� CoronaryInfo ͼ�����ؼ����� CoronaryPara �������ڵ�λ��
% ����ֵ S �������� /tridata �㡢�桢FFR /Validflag �ж������Ƿ�������
[S,tridata,Validflag] = SkeletonGeometry(CoronaryVol,CoronaryInfo,CoronaryPara);
xyz =round(S(:,25:27));
% plot3(P(:,1),P(:,2),P(:,3),'g.','markersize',8,'linewidth',1)
% hold on 
% scatter3(xyz(:,1),xyz(:,2),xyz(:,3),'r.')
dlmwrite('sketeon_test_P.xyz',P,' ');
dlmwrite('sketeon_test.xyz',S(:,25:27),' ');