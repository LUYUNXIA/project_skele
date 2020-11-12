function [S,tridata,Validflag] = SkeletonGeometry(CoronaryVol,CoronaryInfo,CoronaryPara)
 coder.inline('never');
%SkeletonGeometry  generate the skeleton from coronary
%
%   Input arguments:
%       CoronaryVol           - the 3D mask volume for coronary
%       CoronaryInfo          - the info data
%       CoronaryPara          - the parameters for left and right coronary
%                             - a 2XN array
%   Output arguments:
%       S                     - the skeleton array
%   Example:
%       S = SkeletonGeometry(CoronaryVol,CoronaryInfo,CoronaryPara)
%
%
%   Copyright SHENSHI.Inc(杭州晟视科技)
%   Written by H.P. Wang, Date: 2018/11/4
%   Version: 1.0
%   /////////////////////////////////////////////////////////////////
% A的数据结构,包含19列数据
% 龙骨坐标（x y z）   3列
% 龙骨切线（x' y' z'） 3列
% 龙骨二阶导数 (x'', y'', z'') 3列
% 曲率 1列,10
% 曲率半径 1列,11
% 正交面面积 1列,12
% 正交截面近似半径 1列,13
% 血管收缩角,14
% 是否是出口,15
% 父节点的索引,16
% 每个节点的流量,17
% 每个节点的相对压力,18
% 每个节点的FFR,19
% 左右冠脉标志,20
% 是否是分叉,21
% 修改后的直径,22
% 阻抗,23
% 多孔介质区域,24
% 体素坐标下的x,y,z,[25,26,27]
% 体素坐标下的方向,[28,29,30]
% 狭窄率,31
% 备用,[32,33,34,35]
% 出口编号，33
% 未拟合龙骨[36,37,38]
%% resolve the parameters
%入口位置
pos_in = CoronaryPara(:,1:3); 
flag = true(size(pos_in,1),1);
% 对入口坐标进行检查，查找重复的点
for ii = 1:size(pos_in,1)-1
    pos1 = pos_in(ii,1:3);
    for jj = ii+1:size(pos_in,1)
        pos2 = pos_in(jj,1:3);
        dist = sqrt(sum((pos2-pos1).^2,2));
        if dist < 1
            % 重复点
            flag(jj,1) = 0;
        end
    end
end
pos_in = pos_in(flag,1:3);

%入口点的总流量
Q0 = CoronaryPara(1,4);  
%入口压力
P0 = CoronaryPara(:,5)*133;                
w = size(CoronaryVol,1);
l = size(CoronaryVol,2);
h = size(CoronaryVol,3);
Validflag = 0;
% 是否处理末端
CropEnds = CoronaryInfo.CropEnds;

%% 对入口区域进行处理
if size(pos_in,1) == 2
    Vect = pos_in(2,:)-pos_in(1,:);
    Vect = Vect/sqrt(Vect*Vect');
    thre = 20;
    % 注意第二维是xmesh
    % 对入口进行处理
    cur_dir = Vect;
    for ii = 1:2
        cur_pos = pos_in(ii,:);
        if ii == 1
            cur_dir = Vect;
        elseif ii == 2
            cur_dir = -Vect;
        end
        % 沿方向构造中心点
        Cc = cur_pos+thre*cur_dir;
        % 构造矩形区域
        min2 = max(1,round(Cc(1)-0.9*thre));
        max2 = min(l,round(Cc(1)+0.9*thre));
        min1 = max(1,round(Cc(2)-0.9*thre));
        max1 = min(w,round(Cc(2)+0.9*thre));
        min3 = max(1,round(Cc(3)-0.9*thre));
        max3 = min(h,round(Cc(3)+0.9*thre));
        % 将球形区域置0
        [mesh2,mesh1,mesh3] = meshgrid(min2:max2,min1:max1,min3:max3);
        dist = sqrt((mesh1-Cc(2)).^2+(mesh2-Cc(1)).^2+(mesh3-Cc(3)).^2);
        idx0 = find(dist<0.9*thre);
        idx1 = sub2ind([w,l,h],mesh1(idx0),mesh2(idx0),mesh3(idx0));
        CoronaryVol(idx1) = false;
    end
end

%% find the skeleton
skel = Skeleton3D(CoronaryVol);

% 第二维作为x,体素坐标
[skely,skelx,skelz]=ind2sub([w,l,h],find(skel(:)));
skel_num = length(skelx);
A = zeros(skel_num,50);
A(:,1) = skelx;
A(:,2) = skely;
A(:,3) = skelz;
% 控制龙骨点的个数,不能超过10000个
if skel_num>1e4
    % 龙骨生成时间太长，极大可能龙骨有问题
    nodenum = 1;
    S = zeros(nodenum,50);
    tridata.vertex.x = zeros(nodenum,1);
    tridata.vertex.y = zeros(nodenum,1);
    tridata.vertex.z = zeros(nodenum,1);
    tridata.vertex.FFR = zeros(nodenum,1);
    tridata.vertex.StenosisRate = zeros(nodenum,1);
    tridata.faces = zeros(nodenum,3);
    Validflag = -1;
    return;
end


%% 计算龙骨点切向方向
A = gettangentvectors(A);

%% 移除重复点
A = removerepeatednodes(A);

%% 剔除方向不连续的点（孤立点）
nodenum1 = size(A,1);
nodenum2 = 0;
while(nodenum1-nodenum2>0)
    nodenum1 = size(A,1);
    A = removesolitarypoints(A);
    nodenum2 = size(A,1);
end

%% 对点的方向进行平滑
A = smoothdirection(A,1);

%% 根据入口位置定位出入口点,从龙骨点找出离入口最近的点
inlet_num = size(pos_in,1);
for ii = 1:inlet_num
    dis = sqrt((A(:,1)-pos_in(ii,1)).^2+(A(:,2)-pos_in(ii,2)).^2+(A(:,3)-pos_in(ii,3)).^2);
    % 找距离入口最近的点
    [val,ind] = min(dis);
    % 2代表入口，1代表出口
    A(ind,15) = 2;
    % 加入总流量，还需要分配
    A(ind,17) = Q0;
    % 入口压力
    A(ind,18) = P0(ii);
end

%% 对出口进行检查,是否处理末端
A = outletsvalidation(A,CropEnds);

%% 评估龙骨质量
nodenum = size(A,1);
outlet_nodes = find(A(:,15)==1);
routenum = length(outlet_nodes);
if nodenum>4000 || routenum>40
    % 龙骨有问题
    S = A;
    tridata.vertex.x = zeros(nodenum,1);
    tridata.vertex.y = zeros(nodenum,1);
    tridata.vertex.z = zeros(nodenum,1);
    tridata.vertex.FFR = zeros(nodenum,1);
    tridata.vertex.StenosisRate = zeros(nodenum,1);
    tridata.faces = zeros(nodenum,3);
    Validflag = -1;
    return;
end

%% 再次搜索出口到入口，剔除不连通的点，自动建立索引
% 首先清除索引
nodenum1 = size(A,1);
nodenum2 = 0;
while(nodenum1-nodenum2>0)
    nodenum1 = size(A,1);
    A(:,16) = 0;
    [A,isValid] = findconnectedpath(A);
    if isValid == -1
        % 龙骨有问题
        nodenum = size(A,1);
        S = A;
        tridata.vertex.x = zeros(nodenum,1);
        tridata.vertex.y = zeros(nodenum,1);
        tridata.vertex.z = zeros(nodenum,1);
        tridata.vertex.FFR = zeros(nodenum,1);
        tridata.vertex.StenosisRate = zeros(nodenum,1);
        tridata.faces = zeros(nodenum,3);
        Validflag = -1;
        return;
    end   
    nodenum2 = size(A,1);
end

%% 剔除短分支
A = eliminate_small_branches(A);

%% 保存未拟合龙骨坐标
A(:,36:38) = A(:,1:3);

%% 数据拟合,可能会剔除坏点
A = skeletonfitting(A);

%% 数据拟合后有重复点，再次移除重复点
A = removerepeatednodes(A);

%% 重新建立索引关系
nodenum1 = size(A,1);
nodenum2 = 0;
while(nodenum1-nodenum2>0)
    nodenum1 = size(A,1);
    A(:,16) = 0;
    [A,isValid] = findconnectedpath(A);
    if isValid == -1
        % 龙骨有问题
        nodenum = size(A,1);
        S = A;
        tridata.vertex.x = zeros(nodenum,1);
        tridata.vertex.y = zeros(nodenum,1);
        tridata.vertex.z = zeros(nodenum,1);
        tridata.vertex.FFR = zeros(nodenum,1);
        tridata.vertex.StenosisRate = zeros(nodenum,1);
        tridata.faces = zeros(nodenum,3);
        Validflag = -1;
        return;
    end  
    nodenum2 = size(A,1);
end

%% 再次评估龙骨质量
nodenum = size(A,1);
outlet_nodes = find(A(:,15)==1);
routenum = length(outlet_nodes);
if nodenum>3500 || routenum>35
    % 龙骨有问题
    S = A;
    tridata.vertex.x = zeros(nodenum,1);
    tridata.vertex.y = zeros(nodenum,1);
    tridata.vertex.z = zeros(nodenum,1);
    tridata.vertex.FFR = zeros(nodenum,1);
    tridata.vertex.StenosisRate = zeros(nodenum,1);
    tridata.faces = zeros(nodenum,3);
    Validflag = -1;
    return;
end

%% 根据龙骨点匹配冠脉点云
CoronaryVol = matchPCs2Skel(CoronaryVol,A);  

%% 保存对应的体素坐标
A(:,25:30) = A(:,1:6);

%% 计算等值面
% 注意：FV和Nvec中都是第二维为x,以后需要根据这个改变程序
kernel = ones(3,3,3)/27;
CoronaryVol_S = zeros(size(CoronaryVol),'uint8');
CoronaryVol_S(CoronaryVol) = 255;
CoronaryVol_S = imfilter(CoronaryVol_S,kernel,'same','replicate');
% FS,VS只用于显示
[F,V] = MarchingCubes(CoronaryVol_S,255*0.1);
FS.vertices = V;
FS.faces = F;
[F,V] = MarchingCubes(CoronaryVol,0);
FV.vertices = V;
FV.faces = F;

%% 坐标转换
A = convert2physics(A,CoronaryInfo);
FV.vertices = main_for_axis2(FV.vertices,CoronaryInfo);
FS.vertices = main_for_axis2(FS.vertices,CoronaryInfo);

%% 计算截面半径
A = calradius(A,FV);

%% 对半径进行拟合
A = radiusfitting(A);

%% 修正分叉半径
B = RepairRadius(A);
A(:,22)=B(:,13);

%% 分配入口流量
IND = find(A(:,15)==2);
D = 2*A(IND,13);
D = D.^1;
D = D/sum(D);
Q = Q0.*D;
A(IND,17) = Q;

%% find the routes from inlet to outlet
nodenum = size(A,1);
outlet_nodes = find(A(:,15)==1);
routenum = length(outlet_nodes);
% routes 是从入口到出口的路线,存储索引
routes = cell(routenum,1);
for ii = 1:routenum
    TMP = zeros(nodenum,1);
    % 第一个点为出口
    counter = 1;
    TMP(counter,1) = outlet_nodes(ii);
    while (1)
        current_node = TMP(counter,1);
        % 找到当前节点对应对应的父节点
        parent_node = A(current_node,16);
        if parent_node == 0
            break;
        else
            counter = counter+1;
            TMP(counter,1) = parent_node;
        end
    end
    TMP = TMP(1:counter,1);
    % TMP是从出口到入口，需要调整顺序
    TMP = flip(TMP,1);
    routes{ii,1} = TMP;
end

%% index the left and right coronry
% 找到入口点
inlet_nodes = find(A(:,15)==2);
inlet_num = length(inlet_nodes);
% 注意pos_in应该人为给定，以确定每个入口的流量
% pos_in第一个点为左冠，第二个点点为右冠
% 每条分支对应的入口位置
pos_in = main_for_axis2(pos_in,CoronaryInfo);
for ii = 1:routenum
    current_route = routes{ii};
    current_inlet = repmat(A(current_route(1),1:3),inlet_num,1);
    % 查找分支位于哪一个入口
    dist = sqrt((current_inlet(:,1)-pos_in(:,1)).^2+(current_inlet(:,2)-pos_in(:,2)).^2+(current_inlet(:,3)-pos_in(:,3)).^2);
    [val,ind] = min(dist);
    % 标记这条分支的左右属性
    A(current_route,20) = ind;
end

%% 返回数据
S = A;
tridata.vertex.x = FS.vertices(:,1);
tridata.vertex.y = FS.vertices(:,2);
tridata.vertex.z = FS.vertices(:,3);
tridata.vertex.FFR = zeros(length(FS.vertices(:,1)),1);
tridata.vertex.StenosisRate = zeros(length(FS.vertices(:,1)),1);
tridata.faces = FS.faces-1;
end


function skel = Skeleton3D(img,spare)
% SKELETON3D Calculate the 3D skeleton of an arbitrary binary volume using parallel medial axis thinning.
%
% skel = SKELETON3D(img) returns the skeleton of the binary volume 'img'
% skel = SKELETON3D(img,mask) preserves foreground voxels in 'mask'
%
% MATLAB vectorized implementation of the algorithm by Lee, Kashyap and Chu
% "Building skeleton models via 3-D medial surface/axis thinning algorithms."
% Computer Vision, Graphics, and Image Processing, 56(6):462?78, 1994.
%
% Inspired by the ITK implementation by Hanno Homann
% http://hdl.handle.net/1926/1292
% and the Fiji/ImageJ plugin by Ignacio Arganda-Carreras
% http://fiji.sc/wiki/index.php/Skeletonize3D
%
% Philip Kollmannsberger (philipk@gmx.net)
%
% For more information, see <a
% href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/43400-skeleton3d')">Skeleton3D</a> at the MATLAB File Exchange.

% pad volume with zeros to avoid edge effects
skel=padarray(img,[1 1 1]);

if(nargin==2)
    spare=padarray(spare,[1 1 1]);
end

% fill lookup table
eulerLUT = FillEulerLUT;

width = size(skel,1);
height = size(skel,2);
depth = size(skel,3);

unchangedBorders = 0;

while( unchangedBorders < 6 )  % loop until no change for all six border types
    unchangedBorders = 0;
    for currentBorder=1:6 % loop over all 6 directions
        cands = zeros(width,height,depth,'int8');
        %变成int8,不然会变成double
        tmp = int8(skel);
        switch currentBorder
            case 4
                mdiff = diff(tmp,1,1);
                x=2:width; % identify border voxels as candidates
                cands(x,:,:)=mdiff;
            case 3
                mdiff = diff(tmp,1,1);
                x=1:width-1;
                cands(x,:,:)=-mdiff;
            case 1  
                mdiff = diff(tmp,1,2);
                y=2:height;
                cands(:,y,:)=mdiff;
            case 2
                mdiff = diff(tmp,1,2);
                y=1:height-1;
                cands(:,y,:)=-mdiff;
            case 6
                mdiff = diff(tmp,1,3);
                z=2:depth;
                cands(:,:,z)=mdiff;
            case 5
                mdiff = diff(tmp,1,3);
                z=1:depth-1;
                cands(:,:,z)=-mdiff;
        end
        
        % if excluded voxels were passed, remove them from candidates
        if(nargin==2)
            cands = cands.*~spare;
        end
        
        % make sure all candidates are indeed foreground voxels
        cands_ind = find(cands(:)==1 & skel(:));
        %cands = intersect(find(cands(:)==1),find(skel));
        
        noChange = true;
        
        if(~isempty(cands_ind))
            % get subscript indices of candidates
            [x,y,z]=ind2sub([width height depth],cands_ind);
            
            % get 26-neighbourhood of candidates in volume
            nhood = logical(pk_get_nh(skel,cands_ind));
            
            % remove all endpoints (exactly one nb) from list
            di1 = find(sum(nhood,2)==2);
            nhood(di1,:)=[];
            cands_ind(di1)=[];
            x(di1)=[];
            y(di1)=[];
            z(di1)=[];
            
            % remove all non-Euler-invariant points from list
            di2 = find(~p_EulerInv(nhood, eulerLUT'));
            nhood(di2,:)=[];
            cands_ind(di2)=[];
            x(di2)=[];
            y(di2)=[];
            z(di2)=[];
            
            % remove all non-simple points from list
            di3 = find(~p_is_simple(nhood));
            nhood(di3,:)=[];
            cands_ind(di3)=[];
            x(di3)=[];
            y(di3)=[];
            z(di3)=[];
            
            
            % if any candidates left: divide into 8 independent subvolumes
            if(~isempty(x))
                x1 = find(mod(x,2));
                x2 = find(~mod(x,2));
                y1 = find(mod(y,2));
                y2 = find(~mod(y,2));
                z1 = find(mod(z,2));
                z2 = find(~mod(z,2));
                
                ilst = struct('l',{0,0,0,0,0,0,0,0});
                ilst = cell(8,1);
                ilst{1} = intersect(x1,intersect(y1,z1));
                ilst{2} = intersect(x2,intersect(y1,z1));
                ilst{3} = intersect(x1,intersect(y2,z1));
                ilst{4} = intersect(x2,intersect(y2,z1));
                ilst{5} = intersect(x1,intersect(y1,z2));
                ilst{6} = intersect(x2,intersect(y1,z2));
                ilst{7} = intersect(x1,intersect(y2,z2));
                ilst{8} = intersect(x2,intersect(y2,z2));
                
                idx = [];
                % do parallel re-checking for all points in each subvolume
                for i = 1:8
                    if(~isempty(ilst{i}))
                        idx = ilst{i};
                        li = sub2ind([width height depth],x(idx),y(idx),z(idx));
                        skel(li)=0; % remove points
                        nh = logical(pk_get_nh(skel,li));
                        di_rc = find(~p_is_simple(nh));
                        if(~isempty(di_rc)) % if topology changed: revert
                            skel(li(di_rc))=1;
                        else
                            noChange = false; % at least one voxel removed
                        end
                    end
                end
            end
        end
        
        if( noChange )
            unchangedBorders = unchangedBorders + 1;
        end
        
    end
end

% get rid of padded zeros
skel = skel(2:end-1,2:end-1,2:end-1);
end


function LUT = FillEulerLUT

LUT = zeros(1,255);
LUT(1)  =  1;
LUT(3)  = -1;
LUT(5)  = -1;
LUT(7)  =  1;
LUT(9)  = -3;
LUT(11) = -1;
LUT(13) = -1;
LUT(15) =  1;
LUT(17) = -1;
LUT(19) =  1;
LUT(21) =  1;
LUT(23) = -1;
LUT(25) =  3;
LUT(27) =  1;
LUT(29) =  1;
LUT(31) = -1;
LUT(33) = -3;
LUT(35) = -1;
LUT(37) =  3;
LUT(39) =  1;
LUT(41) =  1;
LUT(43) = -1;
LUT(45) =  3;
LUT(47) =  1;
LUT(49) = -1;
LUT(51) =  1;

LUT(53) =  1;
LUT(55) = -1;
LUT(57) =  3;
LUT(59) =  1;
LUT(61) =  1;
LUT(63) = -1;
LUT(65) = -3;
LUT(67) =  3;
LUT(69) = -1;
LUT(71) =  1;
LUT(73) =  1;
LUT(75) =  3;
LUT(77) = -1;
LUT(79) =  1;
LUT(81) = -1;
LUT(83) =  1;
LUT(85) =  1;
LUT(87) = -1;
LUT(89) =  3;
LUT(91) =  1;
LUT(93) =  1;
LUT(95) = -1;
LUT(97) =  1;
LUT(99) =  3;
LUT(101) =  3;
LUT(103) =  1;

LUT(105) =  5;
LUT(107) =  3;
LUT(109) =  3;
LUT(111) =  1;
LUT(113) = -1;
LUT(115) =  1;
LUT(117) =  1;
LUT(119) = -1;
LUT(121) =  3;
LUT(123) =  1;
LUT(125) =  1;
LUT(127) = -1;
LUT(129) = -7;
LUT(131) = -1;
LUT(133) = -1;
LUT(135) =  1;
LUT(137) = -3;
LUT(139) = -1;
LUT(141) = -1;
LUT(143) =  1;
LUT(145) = -1;
LUT(147) =  1;
LUT(149) =  1;
LUT(151) = -1;
LUT(153) =  3;
LUT(155) =  1;

LUT(157) =  1;
LUT(159) = -1;
LUT(161) = -3;
LUT(163) = -1;
LUT(165) =  3;
LUT(167) =  1;
LUT(169) =  1;
LUT(171) = -1;
LUT(173) =  3;
LUT(175) =  1;
LUT(177) = -1;
LUT(179) =  1;
LUT(181) =  1;
LUT(183) = -1;
LUT(185) =  3;
LUT(187) =  1;
LUT(189) =  1;
LUT(191) = -1;
LUT(193) = -3;
LUT(195) =  3;
LUT(197) = -1;
LUT(199) =  1;
LUT(201) =  1;
LUT(203) =  3;
LUT(205) = -1;
LUT(207) =  1;

LUT(209) = -1;
LUT(211) =  1;
LUT(213) =  1;
LUT(215) = -1;
LUT(217) =  3;
LUT(219) =  1;
LUT(221) =  1;
LUT(223) = -1;
LUT(225) =  1;
LUT(227) =  3;
LUT(229) =  3;
LUT(231) =  1;
LUT(233) =  5;
LUT(235) =  3;
LUT(237) =  3;
LUT(239) =  1;
LUT(241) = -1;
LUT(243) =  1;
LUT(245) =  1;
LUT(247) = -1;
LUT(249) =  3;
LUT(251) =  1;
LUT(253) =  1;
LUT(255) = -1;
end


function EulerInv =  p_EulerInv(img,LUT)

% Calculate Euler characteristic for each octant and sum up
eulerChar = zeros(size(img,1),1);
% Octant SWU
n = ones(size(img,1),1);
n(img(:,25)==1) = bitor(n(img(:,25)==1),128);
n(img(:,26)==1) = bitor(n(img(:,26)==1),64);
n(img(:,16)==1) = bitor(n(img(:,16)==1),32);
n(img(:,17)==1) = bitor(n(img(:,17)==1),16);
n(img(:,22)==1) = bitor(n(img(:,22)==1),8);
n(img(:,23)==1) = bitor(n(img(:,23)==1),4);
n(img(:,13)==1) = bitor(n(img(:,13)==1),2);
eulerChar = eulerChar + LUT(n);
% Octant SEU
n = ones(size(img,1),1);
n(img(:,27)==1) = bitor(n(img(:,27)==1),128);
n(img(:,24)==1) = bitor(n(img(:,24)==1),64);
n(img(:,18)==1) = bitor(n(img(:,18)==1),32);
n(img(:,15)==1) = bitor(n(img(:,15)==1),16);
n(img(:,26)==1) = bitor(n(img(:,26)==1),8);
n(img(:,23)==1) = bitor(n(img(:,23)==1),4);
n(img(:,17)==1) = bitor(n(img(:,17)==1),2);
eulerChar = eulerChar + LUT(n);
% Octant NWU
n = ones(size(img,1),1);
n(img(:,19)==1) = bitor(n(img(:,19)==1),128);
n(img(:,22)==1) = bitor(n(img(:,22)==1),64);
n(img(:,10)==1) = bitor(n(img(:,10)==1),32);
n(img(:,13)==1) = bitor(n(img(:,13)==1),16);
n(img(:,20)==1) = bitor(n(img(:,20)==1),8);
n(img(:,23)==1) = bitor(n(img(:,23)==1),4);
n(img(:,11)==1) = bitor(n(img(:,11)==1),2);
eulerChar = eulerChar + LUT(n);
% Octant NEU
n = ones(size(img,1),1);
n(img(:,21)==1) = bitor(n(img(:,21)==1),128);
n(img(:,24)==1) = bitor(n(img(:,24)==1),64);
n(img(:,20)==1) = bitor(n(img(:,20)==1),32);
n(img(:,23)==1) = bitor(n(img(:,23)==1),16);
n(img(:,12)==1) = bitor(n(img(:,12)==1),8);
n(img(:,15)==1) = bitor(n(img(:,15)==1),4);
n(img(:,11)==1) = bitor(n(img(:,11)==1),2);
eulerChar = eulerChar + LUT(n);
% Octant SWB
n = ones(size(img,1),1);
n(img(:,7)==1) = bitor(n(img(:,7)==1),128);
n(img(:,16)==1) = bitor(n(img(:,16)==1),64);
n(img(:,8)==1) = bitor(n(img(:,8)==1),32);
n(img(:,17)==1) = bitor(n(img(:,17)==1),16);
n(img(:,4)==1) = bitor(n(img(:,4)==1),8);
n(img(:,13)==1) = bitor(n(img(:,13)==1),4);
n(img(:,5)==1) = bitor(n(img(:,5)==1),2);
eulerChar = eulerChar + LUT(n);
% Octant SEB
n = ones(size(img,1),1);
n(img(:,9)==1) = bitor(n(img(:,9)==1),128);
n(img(:,8)==1) = bitor(n(img(:,8)==1),64);
n(img(:,18)==1) = bitor(n(img(:,18)==1),32);
n(img(:,17)==1) = bitor(n(img(:,17)==1),16);
n(img(:,6)==1) = bitor(n(img(:,6)==1),8);
n(img(:,5)==1) = bitor(n(img(:,5)==1),4);
n(img(:,15)==1) = bitor(n(img(:,15)==1),2);
eulerChar = eulerChar + LUT(n);
% Octant NWB
n = ones(size(img,1),1);
n(img(:,1)==1) = bitor(n(img(:,1)==1),128);
n(img(:,10)==1) = bitor(n(img(:,10)==1),64);
n(img(:,4)==1) = bitor(n(img(:,4)==1),32);
n(img(:,13)==1) = bitor(n(img(:,13)==1),16);
n(img(:,2)==1) = bitor(n(img(:,2)==1),8);
n(img(:,11)==1) = bitor(n(img(:,11)==1),4);
n(img(:,5)==1) = bitor(n(img(:,5)==1),2);
eulerChar = eulerChar + LUT(n);
% Octant NEB
n = ones(size(img,1),1);
n(img(:,3)==1) = bitor(n(img(:,3)==1),128);
n(img(:,2)==1) = bitor(n(img(:,2)==1),64);
n(img(:,12)==1) = bitor(n(img(:,12)==1),32);
n(img(:,11)==1) = bitor(n(img(:,11)==1),16);
n(img(:,6)==1) = bitor(n(img(:,6)==1),8);
n(img(:,5)==1) = bitor(n(img(:,5)==1),4);
n(img(:,15)==1) = bitor(n(img(:,15)==1),2);
eulerChar = eulerChar + LUT(n);

EulerInv = false(size(eulerChar));
EulerInv(eulerChar==0) = true;
end


function p_is_simple = p_is_simple(N)

% copy neighbors for labeling
n_p = size(N,1);
p_is_simple = ones(1,n_p);

cube = zeros(26,n_p);
cube(1:13,:)=N(:,1:13)';
cube(14:26,:)=N(:,15:27)';

label = 2*ones(1,n_p);

% for all points in the neighborhood
for i=1:26
    
    idx_1 = find(cube(i,:)==1);
    idx_2 = find(p_is_simple);
    idx = intersect(idx_1,idx_2);
    
    if(~isempty(idx))
        
        % start recursion with any octant that contains the point i
        switch( i )
            
            case {1,2,4,5,10,11,13}
                cube(:,idx) = p_oct_label(1, label, cube(:,idx) );
            case {3,6,12,14}
                cube(:,idx) = p_oct_label(2, label, cube(:,idx) );
            case {7,8,15,16}
                cube(:,idx) = p_oct_label(3, label, cube(:,idx) );
            case {9,17}
                cube(:,idx) = p_oct_label(4, label, cube(:,idx) );
            case {18,19,21,22}
                cube(:,idx) = p_oct_label(5, label, cube(:,idx) );
            case {20,23}
                cube(:,idx) = p_oct_label(6, label, cube(:,idx) );
            case {24,25}
                cube(:,idx) = p_oct_label(7, label, cube(:,idx) );
            case 26,
                cube(:,idx) = p_oct_label(8, label, cube(:,idx) );
        end;
        
        label(idx) = label(idx)+1;
        del_idx = find(label>=4);
        
        if(~isempty(del_idx))
            p_is_simple(del_idx) = 0;
        end
    end
end
end



function cube = p_oct_label(octant, label, cube)

% check if there are points in the octant with value 1
if( octant==1 )
    
    % set points in this octant to current label
    % and recurseive labeling of adjacent octants
    idx_1 = find(cube(1,:) == 1);
    if(~isempty(idx_1))
        cube(1,idx_1) = label(idx_1);
    end
    
    idx_2 = find(cube(2,:) == 1);
    if(~isempty(idx_2))
        cube(2,idx_2) = label(idx_2);
        cube(:,idx_2) = p_oct_label(2,label(idx_2),cube(:,idx_2));
    end
    
    idx_4 = find(cube(4,:) == 1);
    if(~isempty(idx_4))
        cube(4,idx_4) = label(idx_4);
        cube(:,idx_4) = p_oct_label(3,label(idx_4),cube(:,idx_4));
    end
    
    idx_5 = find(cube(5,:) == 1);
    if(~isempty(idx_5))
        cube(5,idx_5) = label(idx_5);
        cube(:,idx_5) = p_oct_label(2,label(idx_5),cube(:,idx_5));
        cube(:,idx_5) = p_oct_label(3,label(idx_5),cube(:,idx_5));
        cube(:,idx_5) = p_oct_label(4,label(idx_5),cube(:,idx_5));
    end
    
    idx_10 = find(cube(10,:) == 1);
    if(~isempty(idx_10))
        cube(10,idx_10) = label(idx_10);
        cube(:,idx_10) = p_oct_label(5,label(idx_10),cube(:,idx_10));
    end
    
    idx_11 = find(cube(11,:) == 1);
    if(~isempty(idx_11))
        cube(11,idx_11) = label(idx_11);
        cube(:,idx_11) = p_oct_label(2,label(idx_11),cube(:,idx_11));
        cube(:,idx_11) = p_oct_label(5,label(idx_11),cube(:,idx_11));
        cube(:,idx_11) = p_oct_label(6,label(idx_11),cube(:,idx_11));
    end
    
    idx_13 = find(cube(13,:) == 1);
    if(~isempty(idx_13))
        cube(13,idx_13) = label(idx_13);
        cube(:,idx_13) = p_oct_label(3,label(idx_13),cube(:,idx_13));
        cube(:,idx_13) = p_oct_label(5,label(idx_13),cube(:,idx_13));
        cube(:,idx_13) = p_oct_label(7,label(idx_13),cube(:,idx_13));
    end
    
end

if( octant==2 )
    
    idx_2 = find(cube(2,:) == 1);
    if(~isempty(idx_2))
        cube(2,idx_2) = label(idx_2);
        cube(:,idx_2) = p_oct_label(1,label(idx_2),cube(:,idx_2));
    end
    
    idx_5 = find(cube(5,:) == 1);
    if(~isempty(idx_5))
        cube(5,idx_5) = label(idx_5);
        cube(:,idx_5) = p_oct_label(1,label(idx_5),cube(:,idx_5));
        cube(:,idx_5) = p_oct_label(3,label(idx_5),cube(:,idx_5));
        cube(:,idx_5) = p_oct_label(4,label(idx_5),cube(:,idx_5));
    end
    
    idx_11 = find(cube(11,:) == 1);
    if(~isempty(idx_11))
        cube(11,idx_11) = label(idx_11);
        cube(:,idx_11) = p_oct_label(1,label(idx_11),cube(:,idx_11));
        cube(:,idx_11) = p_oct_label(5,label(idx_11),cube(:,idx_11));
        cube(:,idx_11) = p_oct_label(6,label(idx_11),cube(:,idx_11));
    end
    
    idx_3 = find(cube(3,:) == 1);
    if(~isempty(idx_3))
        cube(3,idx_3) = label(idx_3);
    end
    
    idx_6 = find(cube(6,:) == 1);
    if(~isempty(idx_6))
        cube(6,idx_6) = label(idx_6);
        cube(:,idx_6) = p_oct_label(4,label(idx_6),cube(:,idx_6));
    end
    
    idx_12 = find(cube(12,:) == 1);
    if(~isempty(idx_12))
        cube(12,idx_12) = label(idx_12);
        cube(:,idx_12) = p_oct_label(6,label(idx_12),cube(:,idx_12));
    end
    
    idx_14 = find(cube(14,:) == 1);
    if(~isempty(idx_14))
        cube(14,idx_14) = label(idx_14);
        cube(:,idx_14) = p_oct_label(4,label(idx_14),cube(:,idx_14));
        cube(:,idx_14) = p_oct_label(6,label(idx_14),cube(:,idx_14));
        cube(:,idx_14) = p_oct_label(8,label(idx_14),cube(:,idx_14));
    end
    
end

if( octant==3 )
    
    idx_4 = find(cube(4,:) == 1);
    if(~isempty(idx_4))
        cube(4,idx_4) = label(idx_4);
        cube(:,idx_4) = p_oct_label(1,label(idx_4),cube(:,idx_4));
    end
    
    idx_5 = find(cube(5,:) == 1);
    if(~isempty(idx_5))
        cube(5,idx_5) = label(idx_5);
        cube(:,idx_5) = p_oct_label(1,label(idx_5),cube(:,idx_5));
        cube(:,idx_5) = p_oct_label(2,label(idx_5),cube(:,idx_5));
        cube(:,idx_5) = p_oct_label(4,label(idx_5),cube(:,idx_5));
    end
    
    idx_13 = find(cube(13,:) == 1);
    if(~isempty(idx_13))
        cube(13,idx_13) = label(idx_13);
        cube(:,idx_13) = p_oct_label(1,label(idx_13),cube(:,idx_13));
        cube(:,idx_13) = p_oct_label(5,label(idx_13),cube(:,idx_13));
        cube(:,idx_13) = p_oct_label(7,label(idx_13),cube(:,idx_13));
    end
    
    idx_7 = find(cube(7,:) == 1);
    if(~isempty(idx_7))
        cube(7,idx_7) = label(idx_7);
    end
    
    idx_8 = find(cube(8,:) == 1);
    if(~isempty(idx_8))
        cube(8,idx_8) = label(idx_8);
        cube(:,idx_8) = p_oct_label(4,label(idx_8),cube(:,idx_8));
    end
    
    idx_15 = find(cube(15,:) == 1);
    if(~isempty(idx_15))
        cube(15,idx_15) = label(idx_15);
        cube(:,idx_15) = p_oct_label(7,label(idx_15),cube(:,idx_15));
    end
    
    idx_16 = find(cube(16,:) == 1);
    if(~isempty(idx_13))
        cube(16,idx_16) = label(idx_16);
        cube(:,idx_16) = p_oct_label(4,label(idx_16),cube(:,idx_16));
        cube(:,idx_16) = p_oct_label(7,label(idx_16),cube(:,idx_16));
        cube(:,idx_16) = p_oct_label(8,label(idx_16),cube(:,idx_16));
    end
    
end

if( octant==4 )
    
    idx_5 = find(cube(5,:) == 1);
    if(~isempty(idx_5))
        cube(5,idx_5) = label(idx_5);
        cube(:,idx_5) = p_oct_label(1,label(idx_5),cube(:,idx_5));
        cube(:,idx_5) = p_oct_label(2,label(idx_5),cube(:,idx_5));
        cube(:,idx_5) = p_oct_label(3,label(idx_5),cube(:,idx_5));
    end
    
    idx_6 = find(cube(6,:) == 1);
    if(~isempty(idx_6))
        cube(6,idx_6) = label(idx_6);
        cube(:,idx_6) = p_oct_label(2,label(idx_6),cube(:,idx_6));
    end
    
    idx_14 = find(cube(14,:) == 1);
    if(~isempty(idx_14))
        cube(14,idx_14) = label(idx_14);
        cube(:,idx_14) = p_oct_label(2,label(idx_14),cube(:,idx_14));
        cube(:,idx_14) = p_oct_label(6,label(idx_14),cube(:,idx_14));
        cube(:,idx_14) = p_oct_label(8,label(idx_14),cube(:,idx_14));
    end
    
    idx_8 = find(cube(8,:) == 1);
    if(~isempty(idx_8))
        cube(8,idx_8) = label(idx_8);
        cube(:,idx_8) = p_oct_label(3,label(idx_8),cube(:,idx_8));
    end
    
    idx_16 = find(cube(16,:) == 1);
    if(~isempty(idx_16))
        cube(16,idx_16) = label(idx_16);
        cube(:,idx_16) = p_oct_label(3,label(idx_16),cube(:,idx_16));
        cube(:,idx_16) = p_oct_label(7,label(idx_16),cube(:,idx_16));
        cube(:,idx_16) = p_oct_label(8,label(idx_16),cube(:,idx_16));
    end
    
    idx_9 = find(cube(9,:) == 1);
    if(~isempty(idx_9))
        cube(9,idx_9) = label(idx_9);
    end
    
    idx_17 = find(cube(17,:) == 1);
    if(~isempty(idx_17))
        cube(17,idx_17) = label(idx_17);
        cube(:,idx_17) = p_oct_label(8,label(idx_17),cube(:,idx_17));
    end
    
end

if( octant==5 )
    
    idx_10 = find(cube(10,:) == 1);
    if(~isempty(idx_10))
        cube(10,idx_10) = label(idx_10);
        cube(:,idx_10) = p_oct_label(1,label(idx_10),cube(:,idx_10));
    end
    
    idx_11 = find(cube(11,:) == 1);
    if(~isempty(idx_11))
        cube(11,idx_11) = label(idx_11);
        cube(:,idx_11) = p_oct_label(1,label(idx_11),cube(:,idx_11));
        cube(:,idx_11) = p_oct_label(2,label(idx_11),cube(:,idx_11));
        cube(:,idx_11) = p_oct_label(6,label(idx_11),cube(:,idx_11));
    end
    
    idx_13 = find(cube(13,:) == 1);
    if(~isempty(idx_13))
        cube(13,idx_13) = label(idx_13);
        cube(:,idx_13) = p_oct_label(1,label(idx_13),cube(:,idx_13));
        cube(:,idx_13) = p_oct_label(3,label(idx_13),cube(:,idx_13));
        cube(:,idx_13) = p_oct_label(7,label(idx_13),cube(:,idx_13));
    end
    
    idx_18 = find(cube(18,:) == 1);
    if(~isempty(idx_18))
        cube(18,idx_18) = label(idx_18);
    end
    
    idx_19 = find(cube(19,:) == 1);
    if(~isempty(idx_19))
        cube(19,idx_19) = label(idx_19);
        cube(:,idx_19) = p_oct_label(6,label(idx_19),cube(:,idx_19));
    end
    
    idx_21 = find(cube(21,:) == 1);
    if(~isempty(idx_21))
        cube(21,idx_21) = label(idx_21);
        cube(:,idx_21) = p_oct_label(7,label(idx_21),cube(:,idx_21));
    end
    
    idx_22 = find(cube(22,:) == 1);
    if(~isempty(idx_22))
        cube(22,idx_22) = label(idx_22);
        cube(:,idx_22) = p_oct_label(6,label(idx_22),cube(:,idx_22));
        cube(:,idx_22) = p_oct_label(7,label(idx_22),cube(:,idx_22));
        cube(:,idx_22) = p_oct_label(8,label(idx_22),cube(:,idx_22));
    end
    
end

if( octant==6 )
    
    idx_11 = find(cube(11,:) == 1);
    if(~isempty(idx_11))
        cube(11,idx_11) = label(idx_11);
        cube(:,idx_11) = p_oct_label(1,label(idx_11),cube(:,idx_11));
        cube(:,idx_11) = p_oct_label(2,label(idx_11),cube(:,idx_11));
        cube(:,idx_11) = p_oct_label(5,label(idx_11),cube(:,idx_11));
    end
    
    idx_12 = find(cube(12,:) == 1);
    if(~isempty(idx_12))
        cube(12,idx_12) = label(idx_12);
        cube(:,idx_12) = p_oct_label(2,label(idx_12),cube(:,idx_12));
    end
    
    idx_14 = find(cube(14,:) == 1);
    if(~isempty(idx_14))
        cube(14,idx_14) = label(idx_14);
        cube(:,idx_14) = p_oct_label(2,label(idx_14),cube(:,idx_14));
        cube(:,idx_14) = p_oct_label(4,label(idx_14),cube(:,idx_14));
        cube(:,idx_14) = p_oct_label(8,label(idx_14),cube(:,idx_14));
    end
    
    idx_19 = find(cube(19,:) == 1);
    if(~isempty(idx_19))
        cube(19,idx_19) = label(idx_19);
        cube(:,idx_19) = p_oct_label(5,label(idx_19),cube(:,idx_19));
    end
    
    
    idx_22 = find(cube(22,:) == 1);
    if(~isempty(idx_22))
        cube(22,idx_22) = label(idx_22);
        cube(:,idx_22) = p_oct_label(5,label(idx_22),cube(:,idx_22));
        cube(:,idx_22) = p_oct_label(7,label(idx_22),cube(:,idx_22));
        cube(:,idx_22) = p_oct_label(8,label(idx_22),cube(:,idx_22));
    end
    
    idx_20 = find(cube(20,:) == 1);
    if(~isempty(idx_20))
        cube(20,idx_20) = label(idx_20);
    end
    
    idx_23 = find(cube(23,:) == 1);
    if(~isempty(idx_23))
        cube(23,idx_23) = label(idx_23);
        cube(:,idx_23) = p_oct_label(8,label(idx_23),cube(:,idx_23));
    end
    
end

if( octant==7 )
    
    idx_13 = find(cube(13,:) == 1);
    if(~isempty(idx_13))
        cube(13,idx_13) = label(idx_13);
        cube(:,idx_13) = p_oct_label(1,label(idx_13),cube(:,idx_13));
        cube(:,idx_13) = p_oct_label(3,label(idx_13),cube(:,idx_13));
        cube(:,idx_13) = p_oct_label(5,label(idx_13),cube(:,idx_13));
    end
    
    idx_15 = find(cube(15,:) == 1);
    if(~isempty(idx_15))
        cube(15,idx_15) = label(idx_15);
        cube(:,idx_15) = p_oct_label(3,label(idx_15),cube(:,idx_15));
    end
    
    idx_16 = find(cube(16,:) == 1);
    if(~isempty(idx_16))
        cube(16,idx_16) = label(idx_16);
        cube(:,idx_16) = p_oct_label(3,label(idx_16),cube(:,idx_16));
        cube(:,idx_16) = p_oct_label(4,label(idx_16),cube(:,idx_16));
        cube(:,idx_16) = p_oct_label(8,label(idx_16),cube(:,idx_16));
    end
    
    idx_21 = find(cube(21,:) == 1);
    if(~isempty(idx_21))
        cube(21,idx_21) = label(idx_21);
        cube(:,idx_21) = p_oct_label(5,label(idx_21),cube(:,idx_21));
    end
    
    idx_22 = find(cube(22,:) == 1);
    if(~isempty(idx_22))
        cube(22,idx_22) = label(idx_22);
        cube(:,idx_22) = p_oct_label(5,label(idx_22),cube(:,idx_22));
        cube(:,idx_22) = p_oct_label(6,label(idx_22),cube(:,idx_22));
        cube(:,idx_22) = p_oct_label(8,label(idx_22),cube(:,idx_22));
    end
    
    idx_24 = find(cube(24,:) == 1);
    if(~isempty(idx_24))
        cube(24,idx_24) = label(idx_24);
    end
    
    idx_25 = find(cube(25,:) == 1);
    if(~isempty(idx_25))
        cube(25,idx_25) = label(idx_25);
        cube(:,idx_25) = p_oct_label(8,label(idx_25),cube(:,idx_25));
    end
end

if( octant==8 )
    
    idx_14 = find(cube(14,:) == 1);
    if(~isempty(idx_14))
        cube(14,idx_14) = label(idx_14);
        cube(:,idx_14) = p_oct_label(2,label(idx_14),cube(:,idx_14));
        cube(:,idx_14) = p_oct_label(4,label(idx_14),cube(:,idx_14));
        cube(:,idx_14) = p_oct_label(6,label(idx_14),cube(:,idx_14));
    end
    
    idx_16 = find(cube(16,:) == 1);
    if(~isempty(idx_16))
        cube(16,idx_16) = label(idx_16);
        cube(:,idx_16) = p_oct_label(3,label(idx_16),cube(:,idx_16));
        cube(:,idx_16) = p_oct_label(4,label(idx_16),cube(:,idx_16));
        cube(:,idx_16) = p_oct_label(7,label(idx_16),cube(:,idx_16));
    end
    
    idx_17 = find(cube(17,:) == 1);
    if(~isempty(idx_17))
        cube(17,idx_17) = label(idx_17);
        cube(:,idx_17) = p_oct_label(4,label(idx_17),cube(:,idx_17));
    end
    
    idx_22 = find(cube(22,:) == 1);
    if(~isempty(idx_22))
        cube(22,idx_22) = label(idx_22);
        cube(:,idx_22) = p_oct_label(5,label(idx_22),cube(:,idx_22));
        cube(:,idx_22) = p_oct_label(6,label(idx_22),cube(:,idx_22));
        cube(:,idx_22) = p_oct_label(7,label(idx_22),cube(:,idx_22));
    end
    
    idx_17 = find(cube(17,:) == 1);
    if(~isempty(idx_17))
        cube(17,idx_17) = label(idx_17);
        cube(:,idx_17) = p_oct_label(4,label(idx_17),cube(:,idx_17));
    end
    
    idx_23 = find(cube(23,:) == 1);
    if(~isempty(idx_23))
        cube(23,idx_23) = label(idx_23);
        cube(:,idx_23) = p_oct_label(6,label(idx_23),cube(:,idx_23));
    end
    
    idx_25 = find(cube(25,:) == 1);
    if(~isempty(idx_25))
        cube(25,idx_25) = label(idx_25);
        cube(:,idx_25) = p_oct_label(7,label(idx_25),cube(:,idx_25));
    end
    
    idx_26 = find(cube(26,:) == 1);
    if(~isempty(idx_26))
        cube(26,idx_26) = label(idx_26);
    end
end
end


function nhood = pk_get_nh(img,i)

width = size(img,1);
height = size(img,2);
depth = size(img,3);

[x,y,z]=ind2sub([width height depth],i);

nhood = false(length(i),27);

for xx=1:3
    for yy=1:3
        for zz=1:3
            w=sub2ind([3 3 3],xx,yy,zz);
            idx = sub2ind([width height depth],x+xx-2,y+yy-2,z+zz-2);
            idx = reshape(idx,[],1);
            nhood(:,w)=img(idx);
        end
    end
end
end


function A = removerepeatednodes(A)
%REMOVEREPEATEDNOES remove the repeated nodes in vessel nodes

%   Copyright
%   $Revision: 1.0 $    $Date: 2018/05/19 $

% 保证节点不重复
nodenum = size(A,1);
thre = 0.2;
IND = true(nodenum,1);
for ii = 1:nodenum-1
    current_pos = A(ii,1:3);
    counter = 1;
    for jj = ii+1:nodenum
        tmp = sqrt((A(jj,1)-current_pos(1))^2+(A(jj,2)-current_pos(2))^2+(A(jj,3)-current_pos(3))^2);
        if tmp < thre
            counter = counter+1;
            IND(jj,1) = false;
        end
    end
end
A = A(IND,:);
end


function A = gettangentvectors(A)
% PCA to get the tangent vectors
num = size(A,1);
% 用于计算方向的点的个数
neednum = 7;
if num < neednum
    neednum = num-1;
end
for ii = 1:num
    cur_p = A(ii,1:3);
    dist0 = sqrt((A(:,1)-cur_p(1)).^2+(A(:,2)-cur_p(2)).^2+(A(:,3)-cur_p(3)).^2);
    % 从中找出距离最小的7个点
    [dist1,IND1] = sort(dist0,'ascend');
    IND = IND1(1:neednum);
    pca_p = A(IND,1:3);
    coeff = pca(pca_p);
    tdir = coeff(:,1);
    tdir = tdir/norm(tdir,2);
    A(ii,4:6) = tdir;
end
end


function A = smoothdirection(A, flag)
%REMOVEOUTLIERS remove the outliers of solitary points

%   Copyright
%   $Revision: 1.0 $    $Date: 2018/05/19 $
nodenum = size(A,1);
thre_r = 4;
switch flag
    case 1
        for ii = 1:nodenum
            % 不计算入口节点
            cur_pos = A(ii,1:3);
            % 计算所有点到当前点的距离
            r = sqrt((A(:,1)-cur_pos(1)).^2+(A(:,2)-cur_pos(2)).^2+(A(:,3)-cur_pos(3)).^2);
            % 计算方向间的投影
            tmp = A(r<=thre_r,4:6);
            proj = tmp*(A(ii,4:6))';
            tmp = repmat(sign(proj),1,3).*tmp;
            % 平均值作为当前的方向
            tmp = mean(tmp,1);
            tmp = tmp/sqrt(tmp*tmp');
            A(ii,4:6) = tmp;
        end
    case 2
        for ii = 1:nodenum
            % 不计算入口节点
            cur_pos = A(ii,1:3);
            % 计算所有点到当前点的距离
            r = sqrt((A(:,1)-cur_pos(1)).^2+(A(:,2)-cur_pos(2)).^2+(A(:,3)-cur_pos(3)).^2);
            % 计算方向间的投影
            tmp = A(r<=thre_r,4:6);
            % 平均值作为当前的方向
            tmp = mean(tmp,1);
            tmp = tmp/sqrt(tmp*tmp');
            A(ii,4:6) = tmp;
        end
end
end


function A = findoutletsandinlets(A)
%FINDOUTLETSANDINLETS find the outlets and inlets of vessels

%   Copyright
%   $Revision: 1.0 $    $Date: 2018/05/19 $
nodenum = size(A,1);
thre_r = 4;
thre_proj = sqrt(2)/2;
% 先正反向构造圆，剔除正反向都没有点的点
for ii = 1:nodenum
    % 不计算入口节点和坏点
    if A(ii,15) == 2
        continue;
    end
    cur_pos = A(ii,1:3);
    cur_dir = A(ii,4:6);
    % 对方向进行归一化
    cur_dir = cur_dir/sqrt(cur_dir*cur_dir');
    % 沿方向构造圆,不知道是入口还是出口，需要再反向构造一次
    Cc = cur_pos+thre_r*cur_dir;
    % 计算所有点到圆心的距离
    r = sqrt((A(:,1)-Cc(1)).^2+(A(:,2)-Cc(2)).^2+(A(:,3)-Cc(3)).^2);
    % 判断在这个圆中是否存在点,注意排除自身（0.99）
    IND = find( r<0.99*thre_r, 1);
    counter = 0;
    if isempty(IND)
        % 在方向范围内没有点，为出口或入口
        counter = counter+1;
    end
    % 反向构造一次
    % 沿方向构造圆,不知道是入口还是出口，需要再反向构造一次
    Cc = cur_pos-thre_r*cur_dir;
    % 计算所有点到圆心的距离
    r = sqrt((A(:,1)-Cc(1)).^2+(A(:,2)-Cc(2)).^2+(A(:,3)-Cc(3)).^2);
    % 判断在这个圆中是否存在点,注意排除自身
    IND = find( r<0.99*thre_r, 1);
    if isempty(IND)
        counter = counter+1;
    end
    % 如果counter = 2,说明该点是个单独的点，需要剔除
    if counter == 0
        % 中间点
    elseif counter == 1
        A(ii,15) = 1;
    elseif counter == 2
        A(ii,15) = -1;
    end
end
A(A(:,15)==-1,:) = [];
end


function A = outletsvalidation(A,CropEnds)
%% 剔除前n个末端出口点
if CropEnds == 1
    % 处理末端
    for ii = 1:12
        A(A(:,15)==1,:) = [];
        A = findoutletsandinlets(A);
    end
elseif CropEnds == 0
    % 不处理末端
    for ii = 1:1
        A(A(:,15)==1,:) = [];
        A = findoutletsandinlets(A);
    end
end
outlet_nodes = find(A(:,15)==1);
outpos = A(A(:,15)==1,1:3);
%% 排除重复的出口点
% 相邻出口点的最大距离
thre = 2;
for ii = 1:size(outpos,1)-1
    outpos_1 = outpos(ii,:);
    for jj = ii:size(outpos,1)
        outpos_2 = outpos(jj,:);
        dist = outpos_1-outpos_2;
        dist = sqrt(dist*dist');
        if dist < thre
            % 这两个出口点靠的太近，都剔除
            A(outlet_nodes(ii),15) = -1;
            A(outlet_nodes(jj),15) = -1;
        end
    end
end
% 重新计算出口
A(A(:,15)==-1,:) = [];
A = findoutletsandinlets(A);
%% 对出口进行检查,若能连续找到5个出口则认为该出口是正确的
% 必须迭代两次，一次可能还是有错误的出口
for iter = 1:2
    B = A;
    Lb = 5;
    % 记录原始出口是不是坏点，1表示正确，0表示错误，初始为全部错误
    outlet_nodes = find(A(:,15)==1);
    outlierflag = zeros(size(outlet_nodes,1),Lb+1);
    outlierflag(:,1) = 1;
    outlierpos = zeros(size(outlet_nodes,1),3,Lb+1);
    outlierpos(:,:,1) = A(outlet_nodes,1:3);
    % 记录原始出口位置
    outpos_1 = A(A(:,15)==1,1:3);
    for ii = 1:Lb
        % 首先将出口点去掉
        B(B(:,15)==1,:) = [];
        % 重新计算出口,矩阵大小不变
        B = findoutletsandinlets(B);
        % 记录新的出口位置
        outpos_2 = B(B(:,15)==1,1:3);
        % 在新的出口位置中查找和原始出口位置邻近的点
        for jj = 1:size(outpos_2,1)
            dist = sqrt((outpos_1(:,1)-outpos_2(jj,1)).^2+(outpos_1(:,2)-outpos_2(jj,2)).^2+(outpos_1(:,3)-outpos_2(jj,3)).^2);
            [val,ind] = min(dist);
            if val<thre*ii
                % 该原始出口点正常
                outlierflag(ind,ii+1) = 1;
                % 记录和原始出口对应的分支出口
                outlierpos(ind,1:3,ii+1) = outpos_2(jj,:);
            end
        end
    end
    outlierflag = sum(outlierflag,2);
    % 找出不符合标准的出口点，并剔除
    outlierpos = outlierpos(outlierflag<Lb+1,:,:);
    for ii = 1:size(outlierpos,1)
        for jj = 1:size(outlierpos,3)
            pos = squeeze(outlierpos(ii,:,jj));
            if sqrt(pos*pos')>100*eps
                dist = sqrt((A(:,1)-pos(1)).^2+(A(:,2)-pos(2)).^2+(A(:,3)-pos(3)).^2);
                [val,ind] = min(dist);
                if val<100*eps
                    A(ind,15) = -1;
                end
            end
        end
    end
    % 重新计算出口
    A(A(:,15)==-1,:) = [];
    A = findoutletsandinlets(A);
end
end


function [A,Validflag] = findconnectedpath(A)
%FINDCONNECTEDPATH 查找出口到入口的连通路径
% 算法流程：
% 1. 从出口开始搜索，记录已经经过的点
% 2. 若碰到走不通的地方（出口、坏点），进行倒车，倒车路径认为是坏点，不再行走
% 3. 倒车过程中再次选择路径，注意不在选择坏点
% 4. 直到达到入口

%   Copyright
%   $Revision: 1.0 $    $Date: 2018/05/19 $

nodenum = size(A,1);
Validflag = 0;
% A的数据结构,包含19列数据
% 龙骨坐标（x y z）   3列
% 龙骨切线（x' y' z'） 3列
% 龙骨二阶导数 (x'', y'', z'') 3列
% 曲率 1列,10
% 曲率半径 1列,11
% 正交面面积 1列,12
% 正交截面近似半径 1列,13
% 血管收缩角,14
% 是否是出口,15
% 父节点的索引,16
% 每个节点的流量,17
% 每个节点的相对压力,18
% 每个节点的FFR,19

%% 对出口节点进行方向判断
outlet_nodes = find(A(:,15) == 1);
routenum = length(outlet_nodes);                %分支个数
spacing = 0;                                    %节点之间的距离
outlierInd = true(nodenum,1);                   %记录坏点
counter = 0;
proj_thre = 0.2;                               % 将proj_thre缩小，更加接近最小距离算法，避免错误方向的影响
for kk = 1:routenum
    % 初始点为出口
    cur_ind1 = outlet_nodes(kk);
    cur_pos = A(cur_ind1,1:3);
    cur_dir = A(cur_ind1,4:6);
    % 与当前点最接近的点
    dis = sqrt((A(:,1)-cur_pos(1)).^2+(A(:,2)-cur_pos(2)).^2+(A(:,3)-cur_pos(3)).^2);
    % 排除自身
    dis(cur_ind1,1) = max(dis(:));
    [val, tmp_ind] = min(dis);
    tmp_dir = A(tmp_ind,1:3)-cur_pos;
    tmp_dir = tmp_dir/sqrt(tmp_dir*tmp_dir');
    proj = cur_dir*tmp_dir';
    % 计算两者的夹角
    if proj > 100*eps
        % 当前方向不变
        counter = counter+1;
        opt_dir = cur_dir;
        A(cur_ind1,4:6) = opt_dir;
        A(cur_ind1,16) = tmp_ind;
        spacing = spacing + val;
    elseif proj < -100*eps
        counter = counter+1;
        opt_dir = -cur_dir;
        A(cur_ind1,4:6) = opt_dir;
        A(cur_ind1,16) = tmp_ind;
        spacing = spacing + val;
    else
        % 两者之间差别较大，该节点有问题
        outlierInd(cur_ind1,1) = false;
    end
end
A = A(outlierInd,:);
spacing = spacing/counter;
%% 从出口节点开始模拟开车过程
thre_r = min(10*spacing,5);                          % 在体素坐标下thre_r设为3
nodenum = size(A,1);
usedInd = false(nodenum,1);                          %记录行驶过的点（选择的点）
inlet_nodes = find(A(:,15) == 2);
inlet_num = length(inlet_nodes);
outlet_nodes = find(A(:,15) == 1);
routenum = length(outlet_nodes);                     %分支个数
for kk = 1:routenum
    outlierInd = false(nodenum,1);                   %记录倒车路径（坏点）
    % 初始点为出口
    ind1 = outlet_nodes(kk);
    pos1 = A(ind1,1:3);
    dir1 = A(ind1,4:6);
    usedInd(ind1,1) = true;
    % 存储当前分支上点的索引
    routeind = zeros(nodenum,1); 
    % 分支上点的个数
    lengthroute = 0;
    % 记录循环次数
    numloop = 0;
    while(1)
         numloop = numloop+1;
         % 避免单个分支出现死循环
         if numloop>10000
             Validflag = -1;
             return;
         end             
        % 看离那个入口点最近
        dist = repmat(A(ind1,1:3),inlet_num,1)-A(inlet_nodes,1:3);
        dist = sqrt(dist(:,1).^2+dist(:,2).^2+dist(:,3).^2);
        [val,ind] = min(dist);
        % 对应的入口点
        inlet_pos = A(inlet_nodes(ind),1:3);
        % 计算当前点到入口点的距离
        tmp = inlet_pos-pos1;
        dir_in = tmp/sqrt(tmp*tmp');
        dist_in = sqrt(tmp*tmp');
        % 在一定范围为找点
        tmp_all = A(:,1:3)-repmat(pos1,nodenum,1);
        % 计算每个节点上的投影
        proj = dir1(1)*tmp_all(:,1)+dir1(2)*tmp_all(:,2)+dir1(3)*tmp_all(:,3);
        % 计算距离
        dist_all = sqrt(tmp_all(:,1).*tmp_all(:,1)+tmp_all(:,2).*tmp_all(:,2)+tmp_all(:,3).*tmp_all(:,3));
        % 投影无量纲
        proj = proj./dist_all;
        % 排除已经循环过的点
        dist_all(routeind(1:lengthroute)) = nan;
        % 排除坏点,不排除坏点(已经循环过的点)可能导致分支丢失
        dist_all(outlierInd) = nan;
        % 查找邻近，方向上的点
        for tt = 1:3
            cur_ind2 = find(dist_all<=5*tt & proj>0);
            if isempty(cur_ind2)
            else
                break;
            end
        end
        if isempty(cur_ind2)
            % 当前点有问题，将当前节点加入坏点行列
            outlierInd(ind1,1) = true;
            % 倒车，回退一个节点
            if lengthroute > 0
                ind1 = routeind(lengthroute);
                pos1 = A(ind1,1:3);
                dir1 = A(ind1,4:6);
                lengthroute = lengthroute-1;
                continue;
            elseif lengthroute == 0
                A(ind1,15) = 0;
                break;
            end
        end
        % 从这些邻近点中选出和入口最近的点
        tmp = A(cur_ind2,1:3)-repmat(inlet_pos,length(cur_ind2),1);
        dist = sqrt(tmp(:,1).*tmp(:,1)+tmp(:,2).*tmp(:,2)+tmp(:,3).*tmp(:,3));
        [val,minind] = min(dist);
        short_pos = A(cur_ind2(minind),1:3);
        % 根据最近点预测方向
        dir_pre = short_pos-pos1;
        dir_pre = dir_pre/sqrt(dir_pre*dir_pre');
        % 对三个方向进行融合，找到最可能方向,本身的方向占主导
        % 可能需要添加当前点是否有父节点？
        dir = 40/100*dir1+30/100*dir_in+30/100*dir_pre;
        % dir为行驶方向
        dir = dir/sqrt(dir*dir');
        
        %
%                 h1 = plot3(pos1(1),pos1(2),pos1(3),'rx','markersize',8,'linewidth',2);
%                 h2 = quiver3(pos1(1),pos1(2),pos1(3),dir(1),dir(2),dir(3),10,'r');
%                  waitforbuttonpress;
%                 delete(h1);
%                 delete(h2);
        %
        
        % 在该方向上查找最近点,计算每个节点上的投影
        proj = dir(1)*tmp_all(:,1)+dir(2)*tmp_all(:,2)+dir(3)*tmp_all(:,3);
        % 对投影归一化，切向量的模为1
        proj = proj./dist_all;
        proj(outlierInd) = nan;
        % 在给定角度内查找点
        tmpind = find(proj>=proj_thre & dist_all<=thre_r);
        if isempty(tmpind)                  % 当前点没有找到可能的父节点，需要倒车
            % 将当前节点加入坏点行列
            outlierInd(ind1,1) = true;
            % 倒车，回退一个节点
            if lengthroute > 0
                ind1 = routeind(lengthroute);
                pos1 = A(ind1,1:3);
                dir1 = A(ind1,4:6);
                lengthroute = lengthroute-1;
            elseif lengthroute == 0
            end
        else
            dist_1 = A(tmpind,1)-pos1(1);
            dist_2 = A(tmpind,2)-pos1(2);
            dist_3 = A(tmpind,3)-pos1(3);
            dist = sqrt(dist_1.^2+dist_2.^2+dist_3.^2);
            [val,ind2] = min(dist);
            % 索引ind2对应的节点就是当前节点的父节点
            ind2 = tmpind(ind2);
            pos2 = A(ind2,1:3);
            dir2 = A(ind2,4:6);
            % ind2是出口，重新开始循环
            if A(ind2,15) == 1
                % 找到出口，倒车，回退一个节点
                if lengthroute > 0
                    outlierInd(ind2,1) = true;
                    outlierInd(ind1,1) = true;
                    ind1 = routeind(lengthroute);
                    pos1 = A(ind1,1:3);
                    dir1 = A(ind1,4:6);
                    lengthroute = lengthroute-1;
                elseif lengthroute == 0
                end
            elseif A(ind2,15) == 2                    %如果当前点是入口，则停止搜索
                proj = dir1*dir2';
                if proj < 0
                    A(ind2,4:6) = -dir2;
                end
                usedInd(ind1,1) = true;
                usedInd(ind2,1) = true;
                % 找到一个正确点，加入
                lengthroute = lengthroute+1;
                routeind(lengthroute) = ind1;
                lengthroute = lengthroute+1;
                routeind(lengthroute) = ind2;
                break;
            else
                % 父节点和预测方向之间的偏差不能太大
                proj = dir2*dir';
                if proj > 100*eps
                    % 父节点的方向正确
                    usedInd(ind1,1) = true;
                    lengthroute = lengthroute+1;
                    routeind(lengthroute) = ind1;
                    % 更新ind
                    ind1 = ind2;
                    pos1 = A(ind2,1:3);
                    dir1 = A(ind2,4:6);
                elseif proj < -100*eps
                    % 父节点的方向反了
                    usedInd(ind1,1) = true;
                    lengthroute = lengthroute+1;
                    routeind(lengthroute) = ind1;
                    A(ind2,4:6) = -dir2;
                    % 更新ind1
                    ind1 = ind2;
                    pos1 = A(ind2,1:3);
                    dir1 = A(ind2,4:6);
                else
                    % 两者之间差别较大，该节点有问题
                    outlierInd(ind1,1) = true;
                    outlierInd(ind2,1) = true;
                    ind1 = routeind(lengthroute);
                    pos1 = A(ind1,1:3);
                    dir1 = A(ind1,4:6);
%                     routeind = routeind(1:end-1);
                    lengthroute = lengthroute-1;
                end
            end
        end
        % 回到了起始出口点，将该起始出口点删除
        if lengthroute == 0
            A(ind1,15) = 0;
            break;
        end
    end
    % 入口到出口的路径要足够长
    if lengthroute>30
        % 修正父节点
        son_ind = routeind(1:lengthroute-1);
        parrent_ind = routeind(2:lengthroute);
        % 找到没有父节点的点
        tmpind = find(A(son_ind,16)==0);
        son_ind = son_ind(tmpind);
        parrent_ind = parrent_ind(tmpind);
        A(son_ind,16) = parrent_ind;
    end
end
% 存在父节点的点
A = A(A(:,16)>0 | A(:,15) == 2,:);
end



function A = removesolitarypoints(A)
%REMOVEOUTLIERS remove the outliers of solitary points

%   Copyright
%   $Revision: 1.0 $    $Date: 2018/05/19 $
nodenum = size(A,1);
thre_r = 3;
% 将thre_proj调小，相当于不要用方向作为判据
% thre_proj = sqrt(1)/5;
% for ii = 1:nodenum
%     % 不计算入口节点
%     if A(ii,15) == 2
%         continue;
%     end
%     cur_pos = A(ii,1:3);
%     % 计算所有点到当前点的距离
%     r = sqrt((A(:,1)-cur_pos(1)).^2+(A(:,2)-cur_pos(2)).^2+(A(:,3)-cur_pos(3)).^2);
%     % 剔除本身
%     r(ii) = 10*thre_r;
%     % 计算方向间的投影
%     proj = A(r<=thre_r,4:6)*(A(ii,4:6))';
%     proj = abs(proj);
%     % 判断投影是否太小,若太小排除该点
%     if median(proj)<thre_proj
%         A(ii,15) = -1;
%     end
%     % 周围点的平均位置
%     %     sta_pos = mean(A(r<=thre_r,1:3),1);
%     %     r = cur_pos-sta_pos;
%     %     r = sqrt(r*r');
%     %     if r>0.6*thre_r
%     %         A(ii,15) = -1;
%     %     end
% end
% A(A(:,15)==-1,:) = [];

% nodenum = size(A,1);
% 先正反向构造圆，剔除正反向都没有点的点
for ii = 1:nodenum
    % 不计算入口节点
    if A(ii,15) == 2 || A(ii,15) == 1
        continue;
    end
    cur_pos = A(ii,1:3);
    cur_dir = A(ii,4:6);
    % 对方向进行归一化
    cur_dir = cur_dir/sqrt(cur_dir*cur_dir');
    % 沿方向构造圆,不知道是入口还是出口，需要再反向构造一次
    Cc = cur_pos+thre_r*cur_dir;
    % 计算所有点到圆心的距离
    r = sqrt((A(:,1)-Cc(1)).^2+(A(:,2)-Cc(2)).^2+(A(:,3)-Cc(3)).^2);
    % 判断在这个圆中是否存在点,注意排除自身（0.99）
    IND1 = find(r<0.99*thre_r);
    L1 = length(IND1);
    % 反向构造一次
    % 沿方向构造圆,不知道是入口还是出口，需要再反向构造一次
    Cc = cur_pos-thre_r*cur_dir;
    % 计算所有点到圆心的距离
    r = sqrt((A(:,1)-Cc(1)).^2+(A(:,2)-Cc(2)).^2+(A(:,3)-Cc(3)).^2);
    % 判断在这个圆中是否存在点,注意排除自身
    IND2 = find(r<0.99*thre_r);
    L2 = length(IND2);
    % 如果当前点周围的点太少，剔除当前点
    if L1+L2<=2
        A(ii,15) = -1;
    end
end
A(A(:,15)==-1,:) = [];
end


function raw = main_for_axis2(raw,info)
% 读取dicom中的头文件,根据以下几个字段进行坐标转换
% Image Position:表示图像的左上角在空间坐标系中的x,y,z坐标,单位是毫米.

% Image Orientation:表示图像坐标与解剖学坐标体系对应坐标的夹角余弦值.
% 比如origin为（1,0,0,0,9912，-0.1322）表示x轴与L夹角0°,与PS夹角90°;
% y轴与L夹角90°,与P夹角7.6°(arccos(0.9912)=7.6°),与S夹角97.6°(arccos(-0.1322)=97.6°).
% 表示患者躺在病床上，床有点前倾（即身体有些前倾），观察者从患者脚处正着观察患者.

% Pixel Spacing:表示图像每一张Dicom影像的x方向和y方向上的像素间距,由两个值构成的数组.

% Spacing Between Slices:表示图像每一张Dicom影像的z方向上的像素间距,可以根据Image
%Position相邻两张图像z方向相减得到.

% 上述说明主要参考以下三个博客:
% https://blog.csdn.net/sunyao_123/article/details/72801429
% https://blog.csdn.net/sunyao_123/article/details/78975816
% https://blog.csdn.net/inter_peng/article/details/52097847?locationNum=9&fps=1

% 转换公式: (变换的图像坐标-图像原点坐标) * 像素间距 + 原点坐标

% Input:
% raw为图像坐标的三维点（imshow()中的x,y与实际x,y相反）
% info为对应的dicom切片的头文件信息
% Output:
% 进行坐标转换后,新的坐标,使用国际单位m.

% Example:
% raw = [
%          313,193,40;
%          313,194,40;
%          ..........
%          270,283,40;
%       ];  (size = 799 *3 )
%
%  p_xyz = img2patient(raw,info)

% PS: 若输入为二维图像,则可先通过
% index = logical(reshape(img_ind,512*512,1));  % 获取为1（不为0）的index
% [x,y] = ind2sub([512,512],find(index==1)); % 获取为1（不为0, find(index>0)）位置的坐标
% raw = [x, y, ones(size(x))* z ];   % z为当前切片位置
%
%%
% len = size(raw,1);
%
% iop = info.ImageOrientationPatient;
% ipp = info.ImagePositionPatient;
sbs = info.SpacingBetweenSlices;
ps = info.PixelSpacing;

% b=raw(:,1);
% raw(:,1)=raw(:,2);
% raw(:,2)=b;


%  iop_x = iop(1:3);
%  iop_y = iop(4:6);

% 旋转后的像素间距 矩阵
%  R = [
%      iop_x(1) * ps(1), iop_y(1) * ps(2), 0;
%      iop_x(2) * ps(1), iop_y(2) * ps(2), 0;
%      iop_x(3) * ps(1), iop_y(3) * ps(2), 0;
%      ];

% 原始物理坐标  矩阵; 变换后的所有点需要加上原始物理坐标
% p =[ipp(1), ipp(2), ipp(3) + (-(z_ind-1)) * sbs];  %若没有z坐标,则切片数量*sbs
% p =[ipp(1), ipp(2), ipp(3)];
% P = repmat(p,[len,1]);
% p_z = abs(raw(:,3)-1) * sbs ;
% raw(:,3) = raw(:,3) + p_z;
raw(:,3) = raw(:,3).* sbs ;
raw(:,1:2) = [raw(:,1).*ps(1),raw(:,2).*ps(2)];
% 图像坐标变换
% trans_rwa = [raw(:,2)-1,raw(:,1)-1,zeros(len,1)];
% p_xyz = trans_rwa * R + P;
% 物理坐标单位转换, mm to m
%  p_xyz = p_xyz./1000;
raw(:,1:3) = raw(:,1:3)./1000;
end

function A = calradius(A,FV)
% 根据前后直径来判断当前直径是否正确
nodenum = size(A,1);
vertexnum = size(FV.vertices,1);
IND = false(vertexnum,1);
r0 = 0.003;     %r0为预估的血管半径
dist_thre = 0.0002;   %点到平面的距离
for ii = 1:nodenum
    pos0 = A(ii,1:3);
    dir0 = A(ii,4:6);
    pos0_mat = repmat(pos0,vertexnum,1);
    dir0_mat = repmat(dir0,vertexnum,1);
    cur_r = r0;
    iter = 1;
    cur_r1 = zeros(10,1);
    while iter < 5
        % 平面与边界曲面相减，从边界曲面中选择与平面的交线
        % 选择1.5*cur_r范围内的顶点
        tmp = FV.vertices-pos0_mat;
        dist2center = sqrt(sum(tmp.^2,2));
        % 顶点到当前平面的距离,投影到法向
        tmp = dir0_mat.*tmp;
        dist2plane = abs(sum(tmp,2));
        % 找出满足条件的顶点
        INDLength = 0;
        count = 1;
        while INDLength < 8
            IND = ((dist2center<1.5*cur_r) & (dist2plane<dist_thre*count));
            INDLength = sum(IND);
            count = count+1;
            if count > 4
                break;
            end
        end
        % 得到与平面相交的边界点
        cur_pos = FV.vertices(IND,:);
        dist = cur_pos-repmat(pos0,size(cur_pos,1),1);
        dist1 = sqrt(sum(dist.^2,2));
        cur_r = mean(dist1);
        cur_r1(iter,1) = cur_r;
        if iter>=2 && cur_r1(iter) == cur_r1(iter-1)
            break;
        end
        iter = iter+1;
    end
    if isnan(cur_r) || cur_r>=2*r0
        cur_r = 0.002;
    end
    A(ii,13) = cur_r;
    A(ii,12) = pi*cur_r*cur_r;
end
end

function A = skeletonfitting(A)
%% find the routes from inlet to outlet
nodenum = size(A,1);
outlet_nodes = find(A(:,15) == 1);
routenum = length(outlet_nodes);
% routes 是从入口到出口的路线,存储索引
routes = cell(routenum,1);
for ii = 1:routenum
    TMP = zeros(nodenum,1);
    % 第一个点为出口
    counter = 1;
    TMP(counter,1) = outlet_nodes(ii);
    while (1)
        current_node = TMP(counter,1);
        % 找到当前节点对应对应的父节点
        parent_node = A(current_node,16);
        if parent_node == 0
            break;
        else
            counter = counter+1;
            TMP(counter,1) = parent_node;
        end
    end
    TMP = TMP(1:counter,1);
    % TMP是从出口到入口，需要调整顺序
    TMP = flip(TMP,1);
    routes{ii,1} = TMP;
end

%% curve fitting
B = zeros(nodenum,11);
IND = zeros(nodenum,11);
for ii = 1:routenum
    s = A(routes{ii,1},1:3);
    n = size(s,1);
    tmp = zeros(n,11);
    tmp(:,1:3) = s;
    % 计算梯度和二阶导数
    pnum = 3;                                 % the number of used points
    for kk = 1:2
        for jj = 1:n
            idx1 = max(1,jj-pnum);
            idx2 = min(n,jj+pnum);
            tx = [idx1:idx2]'-jj;
            p = polyfit(tx,tmp(idx1:idx2,1),2);     % second order Polynomial curve fitting
            tmp(jj,1) = polyval(p,0);
            tmp(jj,4) = -p(2);
            tmp(jj,7) = p(1);
            p = polyfit(tx,tmp(idx1:idx2,2),2);     % second order Polynomial curve fitting
            tmp(jj,2) = polyval(p,0);
            tmp(jj,5) = -p(2);
            tmp(jj,8) = p(1);
            p = polyfit(tx,tmp(idx1:idx2,3),2);     % second order Polynomial curve fitting
            tmp(jj,3) = polyval(p,0);
            tmp(jj,6) = -p(2);
            tmp(jj,9) = p(1);
        end
    end
    % 计算曲率半径curvature=(r'xr'')/(r'^3)
    r1 = cross(tmp(:,4:6),tmp(:,7:9),2);
    r1 = sqrt(sum(r1.^2,2));
    r2 = sqrt(sum(tmp(:,4:6).^2,2));
    tmp(:,10) = r1./r2.^3;
    tmp(:,11) = r2.^3./r1;
    
    B(routes{ii,1},:) = B(routes{ii,1},:)+tmp;
    IND(routes{ii,1},:) = IND(routes{ii,1},:)+ones(size(tmp));
end
A(:,1:11) = B./IND;
A(isnan(A(:,1)),:) = [];
% 对切向量归一化处理
nodenum = size(A,1);
for ii = 1:nodenum
    tmp = A(ii,4:6);
    proj = sqrt(tmp*tmp');
    A(ii,4:6) = A(ii,4:6)/proj;
end
end


function A = convert2physics(A,info)
%% cconvert the postion to physical units
A(:,1:3) = main_for_axis2(A(:,1:3),info);
%% find the routes from inlet to outlet
nodenum = size(A,1);
outlet_nodes = find(A(:,15) == 1);
routenum = length(outlet_nodes);
% routes 是从入口到出口的路线,存储索引
routes = cell(routenum,1);
for ii = 1:routenum
    TMP = zeros(nodenum,1);
    % 第一个点为出口
    counter = 1;
    TMP(counter,1) = outlet_nodes(ii);
    while (1)
        current_node = TMP(counter,1);
        % 找到当前节点对应对应的父节点
        parent_node = A(current_node,16);
        if parent_node == 0
            break;
        else
            counter = counter+1;
            TMP(counter,1) = parent_node;
        end
    end
    TMP = TMP(1:counter,1);
    % TMP是从出口到入口，需要调整顺序
    TMP = flip(TMP,1);
    routes{ii,1} = TMP;
end

%% curve fitting
B = zeros(nodenum,11);
IND = zeros(nodenum,11);
for ii = 1:routenum
    s = A(routes{ii,1},1:3);
    n = size(s,1);
    tmp = zeros(n,11);
    tmp(:,1:3) = s;
    % 计算梯度和二阶导数
    pnum = 3;                                   % the number of used points
    for jj = 1:n
        idx1 = max(1,jj-pnum);
        idx2 = min(n,jj+pnum);
        tx = (idx1:idx2)'-jj;
        p = polyfit(tx,s(idx1:idx2,1),2);     % second order Polynomial curve fitting
        tmp(jj,4) = -p(2);
        tmp(jj,7) = p(1);
        p = polyfit(tx,s(idx1:idx2,2),2);     % second order Polynomial curve fitting
        tmp(jj,5) = -p(2);
        tmp(jj,8) = p(1);
        p = polyfit(tx,s(idx1:idx2,3),2);     % second order Polynomial curve fitting
        tmp(jj,6) = -p(2);
        tmp(jj,9) = p(1);
    end
    % 计算曲率半径curvature=(r'xr'')/(r'^3)
    r1 = cross(tmp(:,4:6),tmp(:,7:9),2);
    r1 = sqrt(sum(r1.^2,2));
    r2 = sqrt(sum(tmp(:,4:6).^2,2));
    tmp(:,10) = r1./r2.^3;
    tmp(:,11) = r2.^3./r1;
    %     r1 = sum(tmp(:,7:9).^2,2).*sum(tmp(:,4:6).^2,2);
    %     r2 = (sum(tmp(:,4:6).*tmp(:,7:9),2)).^2;
    %     r3 = (sum(tmp(:,4:6).^2,2)).^1.5;
    %     tmp(:,10) = sqrt(r1-r2)./r3;
    %     tmp(:,11) = r3./sqrt(r1-r2);
    
    B(routes{ii,1},:) = B(routes{ii,1},:)+tmp;
    IND(routes{ii,1},:) = IND(routes{ii,1},:)+ones(size(tmp));
end
IND(IND==0) = 1;
A(:,1:11) = B./IND;
A(isnan(A(:,1)),:) = [];
% 对切向量归一化处理
nodenum = size(A,1);
for ii = 1:nodenum
    tmp = A(ii,4:6);
    proj = sqrt(tmp*tmp');
    A(ii,4:6) = A(ii,4:6)/proj;
end
end


function A = radiusfitting(A)
%% find the routes from inlet to outlet
nodenum = size(A,1);
outlet_nodes = find(A(:,15) == 1);
routenum = length(outlet_nodes);
% routes 是从入口到出口的路线,存储索引
routes = cell(routenum,1);
for ii = 1:routenum
    TMP = zeros(nodenum,1);
    % 第一个点为出口
    counter = 1;
    TMP(counter,1) = outlet_nodes(ii);
    while (1)
        current_node = TMP(counter,1);
        % 找到当前节点对应对应的父节点
        parent_node = A(current_node,16);
        if parent_node == 0
            break;
        else
            counter = counter+1;
            TMP(counter,1) = parent_node;
        end
    end
    TMP = TMP(1:counter,1);
    % TMP是从出口到入口，需要调整顺序
    TMP = flip(TMP,1);
    routes{ii,1} = TMP;
end      

%% radius fitting
B = zeros(nodenum,1);
IND = zeros(nodenum,1);
for ii = 1:routenum
    rad = A(routes{ii,1},13);
    n = size(rad,1);
    % 修复入口直径,入口5个点都设置成最大半径
    pnum = 5;
    idx1 = 1;
    idx2 = 1+pnum;
    rad(idx1:idx2) = max(rad(idx1:idx2));
    % 修复出口直径,出口5个点都设置成最大半径
    pnum = 5;
    idx1 = n-pnum+1;
    idx2 = n;
    rad(idx1:idx2) = max(rad(idx1:idx2));
    % 2次平滑
    pnum = 3;
    for iter = 1:2
        for jj = 1:n
            idx1 = max(1,jj-pnum);
            idx2 = min(n,jj+pnum);
            tx = (idx1:idx2)'-jj;
            p = polyfit(tx,rad(idx1:idx2),2);     % second order Polynomial curve fitting
            rad(jj,1) = polyval(p,0);
        end
    end
    B(routes{ii,1},1) = B(routes{ii,1},1)+rad;
    IND(routes{ii,1},1) = IND(routes{ii,1},1)+ones(size(rad));
end
A(:,13) = B./IND;
A(isnan(A(:,13)),:) = [];
end


function [F,V,col] = MarchingCubes(c,iso,colors)

% [F,V] = MarchingCubes(X,Y,Z,C,ISO)
% [F,V,COL] = MarchingCubes(X,Y,Z,C,ISO,COLORS)
%
% Use marching cubes algorithm to compute a triangulated mesh of the
% isosurface within the 3D matrix of scalar values C at isosurface value
% ISO. The 3D matrices (X,Y,Z) represent a Cartesian, axis-aligned grid
% specifying the points at which the data C is given. These coordinate
% arrays must be in the format produced by Matlab's meshgrid function.
% Output arguments F and V are the face list and vertex list
% of the resulting triangulated mesh. The orientation of the triangles is
% chosen such that the normals point from the higher values to the lower
% values. Optional arguments COLORS ans COLS can be used to produce
% interpolated mesh face colors. For usage, see Matlab's isosurface.m.
% To avoid Out of Memory errors when matrix C is large, convert matrices
% X,Y,Z and C from doubles (Matlab default) to singles (32-bit floats).
%
% Adapted for Matlab by Peter Hammer in 2011 based on an
% Octave function written by Martin Helm <martin@mhelm.de> in 2009
% http://www.mhelm.de/octave/m/marching_cube.m
%
% Revised 30 September, 2011 to add code by Oliver Woodford for removing
% duplicate vertices.

PlotFlag = 0;               % 1=plot isosurface, 0=do not plot
calc_cols = false;
lindex = 4;

[edgeTable, triTable] = GetTables();

if ((nargin ~= 2 && nargin ~= 3) || (nargout ~= 2 && nargout ~= 3))
    error('wrong number of input and/or output arguments');
end

if (ndims(c) ~= 3)
    error('c must be matrices of dim 3');
end

if (~isscalar(iso))
    error('iso needs to be scalar value');
end

if ( nargin == 3 )
    if (size(colors) ~= size(c))
        error( 'color must be matrix of same size as c');
    end
    calc_cols = true;
    lindex = 5;
end

n = size(c) - 1; % number of cubes along each direction of image
c = single(c);

% for each cube, assign which edges are intersected by the isosurface

cc = zeros(n(1),n(2),n(3),'uint16'); % 3d array of 8-bit vertex codes

for ii=1:8                             % loop thru vertices of all cubes
    if ii == 1
        idx = c(1:n(1), 1:n(2), 1:n(3)) > iso;  % which cubes have vtx ii > iso
        cc(idx) = bitset(cc(idx), ii);     % for those cubes, turn bit ii on
    elseif ii == 2
        idx = c(2:n(1)+1, 1:n(2), 1:n(3)) > iso;  % which cubes have vtx ii > iso
        cc(idx) = bitset(cc(idx), ii);     % for those cubes, turn bit ii on
    elseif ii == 3
        idx = c(2:n(1)+1, 2:n(2)+1, 1:n(3)) > iso;  % which cubes have vtx ii > iso
        cc(idx) = bitset(cc(idx), ii);     % for those cubes, turn bit ii on
    elseif ii == 4
        idx = c(1:n(1), 2:n(2)+1, 1:n(3)) > iso;  % which cubes have vtx ii > iso
        cc(idx) = bitset(cc(idx), ii);     % for those cubes, turn bit ii on
    elseif ii == 5
        idx = c(1:n(1), 1:n(2), 2:n(3)+1) > iso;  % which cubes have vtx ii > iso
        cc(idx) = bitset(cc(idx), ii);     % for those cubes, turn bit ii on
    elseif ii == 6
        idx = c(2:n(1)+1, 1:n(2), 2:n(3)+1) > iso;  % which cubes have vtx ii > iso
        cc(idx) = bitset(cc(idx), ii);     % for those cubes, turn bit ii on
    elseif ii == 7
        idx = c(2:n(1)+1, 2:n(2)+1, 2:n(3)+1) > iso;  % which cubes have vtx ii > iso
        cc(idx) = bitset(cc(idx), ii);     % for those cubes, turn bit ii on
    elseif ii == 8
        idx = c(1:n(1), 2:n(2)+1, 2:n(3)+1) > iso;  % which cubes have vtx ii > iso
        cc(idx) = bitset(cc(idx), ii);     % for those cubes, turn bit ii on
    end
end
cedge = edgeTable(cc+1);  % intersected edges for each cube ([n1 x n2 x n3] mtx)
id =  find(cedge);        % voxels which are intersected (col of indcs into cedge)
if isempty(id)            % all voxels are above or below iso
    F = [];
    V = [];
    col = [];
    return
end


% calculate the list of intersection points
xyz_off = [1, 1, 1; 2, 1, 1; 2, 2, 1; 1, 2, 1; 1, 1, 2;  2, 1, 2; 2, 2, 2; 1, 2, 2];
edges = [1 2; 2 3; 3 4; 4 1; 5 6; 6 7; 7 8; 8 5; 1 5; 2 6; 3 7; 4 8];
offset = sub2ind(size(c), xyz_off(:, 1), xyz_off(:, 2), xyz_off(:, 3)) -1;
pp = zeros(length(id), lindex, 12);
ccedge = [cedge(id), id];
ix_offset=0;
for jj=1:12
    id__ = logical(bitget(ccedge(:, 1), jj)); % used for logical indexing
    id_ = ccedge(id__, 2);
    [ix,iy,iz] = ind2sub(size(cc), id_);
    id_c = sub2ind(size(c), ix, iy, iz);
    id1 = id_c + offset(edges(jj, 1));
    id2 = id_c + offset(edges(jj, 2));
    % 修改,索引就是坐标
    [ty1,tx1,tz1] = ind2sub(size(c),id1);
    [ty2,tx2,tz2] = ind2sub(size(c),id2);
    if ( calc_cols )
        pp(id__, 1:5, jj) = [InterpolateVertices(iso, tx1, ty1, tz1, ...
            tx2, ty2, tz2, c(id1), c(id2), colors(id1), colors(id2)), ...
            (1:size(id_, 1))' + ix_offset ];
    else
        pp(id__, 1:4, jj) = [InterpolateVertices(iso, tx1, ty1, tz1, ...
            tx2, ty2, tz2, c(id1), c(id2)), ...
            (1:size(id_, 1))' + ix_offset ];
    end
    ix_offset = ix_offset + size(id_, 1);
end

% calculate the triangulation from the point list
totalnum = n*n';
F = zeros(totalnum,3,'double');
counter = 1;
tri = triTable(cc(id)+1, :);
for jj=1:3:15
    id_ = find(tri(:, jj)>0);
    V = [id_, lindex*ones(size(id_, 1), 1),tri(id_, jj:jj+2) ];
    if ( ~ isempty(V) )
        p1 = sub2ind(size(pp), V(:,1), V(:,2), V(:,3));
        p2 = sub2ind(size(pp), V(:,1), V(:,2), V(:,4));
        p3 = sub2ind(size(pp), V(:,1), V(:,2), V(:,5));
        tmp = [pp(p1), pp(p2), pp(p3)];
        F(counter:counter+size(tmp,1)-1,:) = tmp;
        counter = counter+size(tmp,1);
    end
end
F(counter:end,:) = [];

V = zeros(totalnum,3,'double');
col = zeros(totalnum,1,'double');
for jj = 1:12
    idp = pp(:, lindex, jj) > 0;
    if any(idp)
        V(pp(idp, lindex, jj), 1:3) = pp(idp, 1:3, jj);
        if (calc_cols)
            col(pp(idp, lindex, jj),1) = pp(idp, 4, jj);
        end
    end
end
tmp = sum(V.^2,2);
V(tmp==0,:) = [];
col(col==0,:) = [];

% Remove duplicate vertices (by Oliver Woodford)
[V,I] = sortrows(V);
M = [true; any(diff(V), 2)];
V = V(M,:);
I(I) = cumsum(M);
F = I(F);

if PlotFlag
    %     figure('color',[1 1 1])
    %     patch('vertices',V,'faces',F,'edgecolor','none',...
    %         'facecolor',[1 0 0],'facelighting','phong')
    %     light
    %     axis equal off
end
end

% ============================================================
% ==================  SUBFUNCTIONS ===========================
% ============================================================

function p = InterpolateVertices(isolevel,p1x,p1y,p1z,p2x,p2y,p2z,valp1,valp2,col1,col2)

if nargin == 9
    p = zeros(length(p1x), 3);
elseif nargin == 11
    p = zeros(length(p1x), 4);
else
    error('Wrong number of arguments');
end
mu = zeros(length(p1x), 1, 'single');
id = abs(valp1-valp2) < (10*eps) .* (abs(valp1) + abs(valp2));
if ( any(id) )
    p(id, 1:3) = [ p1x(id), p1y(id), p1z(id) ];
    if ( nargin == 11 )
        p(id, 4) = col1(id);
    end
end
nid = ~id;
if any(nid)
    mu(nid) = (isolevel - valp1(nid)) ./ (valp2(nid) - valp1(nid));
    p(nid, 1:3) = [p1x(nid) + mu(nid) .* (p2x(nid) - p1x(nid)), ...
        p1y(nid) + mu(nid) .* (p2y(nid) - p1y(nid)), ...
        p1z(nid) + mu(nid) .* (p2z(nid) - p1z(nid))];
    if nargin == 11
        p(nid, 4) = col1(nid) + mu(nid) .* (col2(nid) - col1(nid));
    end
end
end


function [edgeTable, triTable] = GetTables()

edgeTable = [
    0,     265,  515,  778, 1030, 1295, 1541, 1804, ...
    2060, 2309, 2575, 2822, 3082, 3331, 3593, 3840, ...
    400,   153,  915,  666, 1430, 1183, 1941, 1692, ...
    2460, 2197, 2975, 2710, 3482, 3219, 3993, 3728, ...
    560,   825,   51,  314, 1590, 1855, 1077, 1340, ...
    2620, 2869, 2111, 2358, 3642, 3891, 3129, 3376, ...
    928,   681,  419,  170, 1958, 1711, 1445, 1196, ...
    2988, 2725, 2479, 2214, 4010, 3747, 3497, 3232, ...
    1120, 1385, 1635, 1898,  102,  367,  613,  876, ...
    3180, 3429, 3695, 3942, 2154, 2403, 2665, 2912, ...
    1520, 1273, 2035, 1786,  502,  255, 1013,  764, ...
    3580, 3317, 4095, 3830, 2554, 2291, 3065, 2800, ...
    1616, 1881, 1107, 1370,  598,  863,   85,  348, ...
    3676, 3925, 3167, 3414, 2650, 2899, 2137, 2384, ...
    1984, 1737, 1475, 1226,  966,  719,  453,  204, ...
    4044, 3781, 3535, 3270, 3018, 2755, 2505, 2240, ...
    2240, 2505, 2755, 3018, 3270, 3535, 3781, 4044, ...
    204,   453,  719,  966, 1226, 1475, 1737, 1984, ...
    2384, 2137, 2899, 2650, 3414, 3167, 3925, 3676, ...
    348,    85,  863,  598, 1370, 1107, 1881, 1616, ...
    2800, 3065, 2291, 2554, 3830, 4095, 3317, 3580, ...
    764,  1013,  255,  502, 1786, 2035, 1273, 1520, ...
    2912, 2665, 2403, 2154, 3942, 3695, 3429, 3180, ...
    876,   613,  367,  102, 1898, 1635, 1385, 1120, ...
    3232, 3497, 3747, 4010, 2214, 2479, 2725, 2988, ...
    1196, 1445, 1711, 1958,  170,  419,  681,  928, ...
    3376, 3129, 3891, 3642, 2358, 2111, 2869, 2620, ...
    1340, 1077, 1855, 1590,  314,   51,  825,  560, ...
    3728, 3993, 3219, 3482, 2710, 2975, 2197, 2460, ...
    1692, 1941, 1183, 1430,  666,  915,  153,  400, ...
    3840, 3593, 3331, 3082, 2822, 2575, 2309, 2060, ...
    1804, 1541, 1295, 1030,  778,  515,  265,    0];

triTable =[
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1;
    3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1;
    3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1;
    3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1;
    9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1;
    1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1;
    9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1;
    2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1;
    8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1;
    9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1;
    4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1;
    3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1;
    1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1;
    4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1;
    4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1;
    9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1;
    1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1;
    5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1;
    2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1;
    9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1;
    0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1;
    2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1;
    10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1;
    4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1;
    5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1;
    5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1;
    9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1;
    0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1;
    1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1;
    10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1;
    8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1;
    2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1;
    7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1;
    9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1;
    2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1;
    11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1;
    9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1;
    5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1;
    11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1;
    11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1;
    1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1;
    9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1;
    5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1;
    2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1;
    0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1;
    5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1;
    6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1;
    0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1;
    3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1;
    6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1;
    5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1;
    1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1;
    10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1;
    6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1;
    1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1;
    8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1;
    7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1;
    3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1;
    5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1;
    0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1;
    9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1;
    8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1;
    5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1;
    0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1;
    6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1;
    10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1;
    10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1;
    8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1;
    1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1;
    3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1;
    0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1;
    10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1;
    0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1;
    3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1;
    6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1;
    9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1;
    8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1;
    3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1;
    6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1;
    0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1;
    10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1;
    10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1;
    1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1;
    2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1;
    7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1;
    7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1;
    2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1;
    1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1;
    11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1;
    8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1;
    0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1;
    7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1;
    10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1;
    2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1;
    6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1;
    7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1;
    2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1;
    1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1;
    10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1;
    10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1;
    0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1;
    7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1;
    6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1;
    8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1;
    9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1;
    6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1;
    1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1;
    4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1;
    10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1;
    8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1;
    0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1;
    1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1;
    8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1;
    10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1;
    4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1;
    10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1;
    5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1;
    11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1;
    9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1;
    6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1;
    7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1;
    3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1;
    7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1;
    9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1;
    3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1;
    6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1;
    9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1;
    1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1;
    4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1;
    7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1;
    6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1;
    3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1;
    0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1;
    6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1;
    1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1;
    0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1;
    11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1;
    6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1;
    5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1;
    9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1;
    1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1;
    1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1;
    10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1;
    0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1;
    5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1;
    10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1;
    11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1;
    0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1;
    9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1;
    7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1;
    2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1;
    8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1;
    9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1;
    9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1;
    1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1;
    9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1;
    9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1;
    5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1;
    0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1;
    10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1;
    2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1;
    0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1;
    0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1;
    9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1;
    5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1;
    3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1;
    5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1;
    8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1;
    0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1;
    9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1;
    0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1;
    1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1;
    3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1;
    4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1;
    9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1;
    11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1;
    11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1;
    2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1;
    9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1;
    3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1;
    1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1;
    4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1;
    4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1;
    0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1;
    3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1;
    3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1;
    0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1;
    9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1;
    1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 ] + 1;
end


function CoronaryVol = matchPCs2Skel(CoronaryVol,A)
% 注意A中的X为三维矩阵中的第二维
sizD = size(CoronaryVol);
tmp = false(sizD);
nodenum = size(A,1);
% thre 要尽量大一些
thre = 15;
for ii = 1:nodenum
    Cc = A(ii,1:3);
    % 构造矩形区域
    min2 = max(1,round(Cc(1)-thre));
    max2 = min(sizD(2),round(Cc(1)+thre));
    min1 = max(1,round(Cc(2)-thre));
    max1 = min(sizD(1),round(Cc(2)+thre));
    min3 = max(1,round(Cc(3)-thre));
    max3 = min(sizD(3),round(Cc(3)+thre));
    % 将矩形区域置0
    tmp(min1:max1,min2:max2,min3:max3) = CoronaryVol(min1:max1,min2:max2,min3:max3);
end
CoronaryVol = tmp;
% 对出口进行匹配
thre = 20;
% outlet_nodes = find(A(:,15) == 1);
% outlet_num = length(outlet_nodes);
% for ii = 1:outlet_num
%     cur_pos = A(outlet_nodes(ii),1:3);
%     cur_dir = A(outlet_nodes(ii),4:6);
%     Cc = cur_pos-thre*cur_dir;
%     % 构造球形区域,注意反向
%     Radi = sqrt((A(:,1)-Cc(1)).^2+(A(:,2)-Cc(2)).^2+(A(:,3)-Cc(3)).^2);
%     numinRad = sum(Radi<0.8*thre);
%     % 出口区域必须没有龙骨点
%     if numinRad == 0
%         min2 = max(1,round(Cc(1)-0.8*thre));
%         max2 = min(sizD(2),round(Cc(1)+0.8*thre));
%         min1 = max(1,round(Cc(2)-0.8*thre));
%         max1 = min(sizD(1),round(Cc(2)+0.8*thre));
%         min3 = max(1,round(Cc(3)-0.8*thre));
%         max3 = min(sizD(3),round(Cc(3)+0.8*thre));
%         [mesh2,mesh1,mesh3] = meshgrid(min2:max2,min1:max1,min3:max3);
%         dist = sqrt((mesh1-Cc(2)).^2+(mesh2-Cc(1)).^2+(mesh3-Cc(3)).^2);
%         % 将球形区域置0
%         idx0 = find(dist<0.8*thre);
%         idx1 = sub2ind(sizD,mesh1(idx0),mesh2(idx0),mesh3(idx0));
%         CoronaryVol(idx1) = false;
%     end
% end
% 对入口进行匹配
inlet_nodes = find(A(:,15) == 2);
inlet_num = length(inlet_nodes);
for ii = 1:inlet_num
    cur_pos = A(inlet_nodes(ii),1:3);
    cur_dir = A(inlet_nodes(ii),4:6);
    Cc = cur_pos+thre*cur_dir;
    % 构造球形区域,注意反向
    Radi = sqrt((A(:,1)-Cc(1)).^2+(A(:,2)-Cc(2)).^2+(A(:,3)-Cc(3)).^2);
    numinRad = sum(Radi<0.9*thre);
    % 入口区域必须没有龙骨点
    if numinRad == 0
        min2 = max(1,round(Cc(1)-0.9*thre));
        max2 = min(sizD(2),round(Cc(1)+0.9*thre));
        min1 = max(1,round(Cc(2)-0.9*thre));
        max1 = min(sizD(1),round(Cc(2)+0.9*thre));
        min3 = max(1,round(Cc(3)-0.9*thre));
        max3 = min(sizD(3),round(Cc(3)+0.9*thre));
        [mesh2,mesh1,mesh3] = meshgrid(min2:max2,min1:max1,min3:max3);
        dist = sqrt((mesh1-Cc(2)).^2+(mesh2-Cc(1)).^2+(mesh3-Cc(3)).^2);
        % 将球形区域置0
        idx0 = find(dist<0.9*thre);
        idx1 = sub2ind(sizD,mesh1(idx0),mesh2(idx0),mesh3(idx0));
        CoronaryVol(idx1) = false;
    end
end
end



function B = RepairRadius(A)
kk=1;
c = zeros(40,1);
j=1;
% c =[];
B = A;
outlet_nodes1 = find(B(:,15) == 1);
index = B(:,16);
while kk <= length(outlet_nodes1)
    %     son = [];
    son = zeros(1000,1);
    father1=zeros(1000,1);
    %     father1=[];
    son(1)= outlet_nodes1(kk);
    father1(1)=B(son(1),16);
    i =0;
    while length(find(index==father1(i+1)))==1 && B(son(i+1),15) ~= 2
        i = i+1;
        son(i+1)=father1(i);
        father1(i+1)=B(son(i+1),16);
    end
    if length(find(index==father1(i+1)))==2 && father1(i+1) ~= 0  &&   i>=70
        B_r = B(son(1:i+1),13);  %
        B_r = sort(B_r,'ascend');
        b = 1.2*mean(B_r(i-40:i+1));
        B(son(i-2:i+1),13)=b;
        if  father1(i+1) == 0 || ismember(father1(i+1),c(1:j))
            kk = kk+1;
            continue;
        elseif ~ismember(son(i+1),c(1:j)) && father1(i+1)~= 0
            outlet_nodes1(kk) = father1(i+1);
            c(j) = father1(i+1);
            j=j+1;
            kk = kk-1;
        end
    elseif length(find(index==father1(i+1)))==2 && father1(i+1) ~= 0  && i>=8 && i<20
        if abs(B(son(4),13)-B(son(i),13))> 1.02e-04   %存上往下移
            b = 1.05*mean(B(son(1:i+1),13));
            B(son(i-1:i+1),13)=b;
            if  father1(i+1) == 0 || ismember(father1(i+1),c(1:j))
                kk = kk+1;
                continue;
            elseif ~ismember(son(i+1),c(1:j)) && father1(i+1)~= 0
                outlet_nodes1(kk) = father1(i+1);
                c(j) = father1(i+1);
                j=j+1;
                kk = kk-1;
            end
        end
    elseif length(find(index==father1(i+1)))==2 && father1(i+1) ~= 0  && i>=20 && i<70
        b = 1.2*mean(B(son(3:i),13));
        B(son(i-1:i+1),13)=b;
        if  father1(i+1) == 0 || ismember(father1(i+1),c(1:j))
            kk = kk+1;
            continue;
        elseif ~ismember(son(i+1),c(1:j)) && father1(i+1)~= 0
            outlet_nodes1(kk) = father1(i+1);
            c(j) = father1(i+1);
            j=j+1;
            kk = kk-1;
        end
    elseif  length(find(index==father1(i+1)))==2 && father1(i+1) ~= 0  && i<8 % length(son)< 8
        if  father1(i+1) == 0 || ismember(father1(i+1),c(1:j))
            kk = kk+1;
            continue;
        elseif ~ismember(son(i+1),c(1:j)) && father1(i+1)~= 0
            outlet_nodes1(kk) = father1(i+1);
            c(j) = father1(i+1);
            j=j+1;
            kk = kk-1;
        end
    end
    kk = kk+1;
end
end



function A = eliminate_small_branches(A)
outlet_nodes = find(A(:,15) == 1);
index = A(:,16);        % 读取父节点数据
kk=1;
m=0;
son1 = zeros(22*length(outlet_nodes),1);
while kk <= length(outlet_nodes)
    son = zeros(22,1);
    son(1)= outlet_nodes(kk);
    father=A(son(1),16);
    for i = 1:20
        son(i+1)=father(i);
        father(i+1)=A(son(i+1),16);
        if A(son(i+1),15)==2   %入口处的小分支小于20的去掉
            son(i+1)=0;
            son1(22*m+1:22*(m+1))=son;
            m=m+1;
            break;
        elseif (length(find(index==son(i+1)))==2) && (father(i+1) ~= 0)
            if i<5  %从出口遍历到分叉出点数小于9个点去掉
                son(i+1)=0;
                son1(22*m+1:22*(m+1))=son;
                m=m+1;
                break;
            end
        end
    end
    kk = kk+1;
end
son1 = son1(son1~=0);
% 这些点不是正确的出口
A(son1,15)= 0;
end
