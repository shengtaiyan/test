
close all
clc
clear
dbstop error;
tic
addpath(genpath('distmesh'))
addpath(genpath('PolyMesher'))





fd=@(p) max(ddiff(dcircle(p,0,0,5),dcircle(p,0,0,1)),-min(p(:,1),p(:,2)));
[p0,t0]=distmesh2d(fd,@huniform,1,[0,0;5,5],[0,1;0,5;5,0;1,0;0,1])



% fh = @(p) 0.5-0.25*fd(p)
% [p0,t0]=distmesh2d(fd,fh,0.58,[0,0;5,5],[0,1;0,5;5,0;1,0;0,1])


tricoord0x0      =p0(:,1)
tricoord0y0      =p0(:,2)
tricoord00       =[tricoord0x0 tricoord0y0]
trisdcoon0       =[t0(:,1) t0(:,2) t0(:,3)]
figure
plot(tricoord0x0,tricoord0y0,'r*')

for i=1:length(p0)
    str = i
    text(p0(i,1)+0.1,p0(i,2)+0.1,num2str(str))

end
axis equal
title('tri-point')
tx   =   tricoord0x0(trisdcoon0);
ty   =   tricoord0y0(trisdcoon0);
figure
for itri=1:size(trisdcoon0,1)
    x   =   tricoord0x0(trisdcoon0(itri,:));
    y   =   tricoord0y0(trisdcoon0(itri,:));
    patch(x,y,'white')
    axis equal
    hold on
end
title('tri-mesh')
mat.D                       = IsoElasMtrx(1, 0.3);
mat.den                     = 2
p                           = tricoord00
t                           = trisdcoon0
[ coord0,  sdConn0, sdsc0 ] = triToSBFEMeshWithPlot( p, t )
figure
opt        = struct('LineSpec','-k', 'sdsc0',sdsc0);
PlotSBFEMesh(coord0, sdConn0, opt);
for i=1:length(coord0)
     
    str = i
    text(coord0(i,1)+0.1,coord0(i,2)+0.1,num2str(str))
end
title('MESH0');
coord0
sdConn0
sdsc0



coord0(38,:)

theta38=atan(coord0(38,2)/coord0(38,1))
x38=cos(theta38)
y38=sin(theta38)
coord0(38,:)=[x38 y38]


theta43=atan(coord0(43,2)/coord0(43,1))
x43=cos(theta43)
y43=sin(theta43)
coord0(43,:)=[x43 y43]

theat72=0.5*(theta38+theta43)


x72=cos(theat72)
y72=sin(theat72)




coord0(72,:)=[x72 y72]
sdConn0{6}=[ 43    15
             15    34
             34    32
             32    38 
             38    72
             72    43]
figure
PlotSBFEMesh(coord0, sdConn0, opt);
for i=1:size(sdConn0,1)
    sdsc0(i,:)
    sdsc0(i,1)
    sdsc0(i,2)
    str = i
    text(sdsc0(i,1),sdsc0(i,2),num2str(str))
end
text(coord0(61,1)+0.1,coord0(61,2)+0.1,'A')
text(coord0(71,1)+0.1,coord0(71,2)+0.1,'B')

title('MESH-point');
figure
for i=1:length(coord0)

    str = i
    text(coord0(i,1),coord0(i,2),num2str(str))
end

axis([0,6,0,6])
title('point1');



% figure
% for i=1:size(sdConn0,1)
%     sdsc0(i,:)
%     sdsc0(i,1)
%     sdsc0(i,2)
%     str = i
%     text(sdsc0(i,1),sdsc0(i,2),num2str(str))
% end
% axis([0,6,0,6])
% title('MESH-point2');


figure
opt        = struct('LineSpec','-k', 'sdsc0',sdsc0);
PlotSBFEMesh(coord0, sdConn0, opt);
for i=1:length(coord0)
    str = i;
    text(coord0(i,1)+0.1,coord0(i,2)+0.1,num2str(str))
end
axis([0,5,0,5])
title('MESH0');


figure
opt        = struct('LineSpec','-k', 'sdsc0',sdsc0);
PlotSBFEMesh(coord0, sdConn0, opt);
for i=1:length(sdConn0)
    str = i;
    text(sdsc0(i,1)+0.1,sdsc0(i,2)+0.1,num2str(str))
end
axis([0,5,0,5])
title('MESH0');





figure
opt        = struct('LineSpec','-k', 'sdsc0',sdsc0);
PlotSBFEMesh(coord0, sdConn0, opt);
axis([0,5,0,5])
title('MESH0');

% 
% % C = {1, 2, 3, 4, 5};  % 示例 cell 数组
% target = 71          % 要查找的数字
% 
% % 示例 cell 数组，每个元素是一个矩阵
% C = sdConn0
% 
% % 初始化变量
% found = false;
% matrix_index = 0;
% row_index = 0;
% col_index = 0;
% % 遍历 cell 数组中的每个矩阵
% for i = 1:length(C)
%     % 检查当前矩阵中是否存在目标数字
%     [row, col] = find(C{i} == target);
% 
%     % 如果找到目标数字
%     if ~isempty(row)
%         found = true;
%         matrix_index = i;
%         row_index = row(1);  % 取第一个匹配的位置
%         col_index = col(1);  % 取第一个匹配的位置
%         break;
%     end
% end
% % 输出结果
% if found
%     fprintf('数字 %d 找到了！\n', target);
%     fprintf('它位于 cell 数组的第 %d 个矩阵中，位置是 (%d, %d)。\n', ...
%         matrix_index, row_index, col_index);
% else
%     fprintf('数字 %d 没有找到。\n', target);
% end

a0esong1step1
a0esong1step2
a0esong1step3
a0esong1step4
a0esong1step5
a0esong1step6

a0esong1exactsolution
toc
disp(['运行时间：',num2str(toc)])

 