%a0esong1step5
% displacement boudary condition 位移边界条件
eps          =    1d-5
dispbc4x4    =    find(abs(coord4(:,2))<eps );
dispbc4y4    =    find(abs(coord4(:,1))<eps );
BC_Dispx4    =    [dispbc4x4     2*  ones(size(dispbc4x4,1),1)  zeros(size(dispbc4x4,1),1)];
BC_Dispy4    =    [dispbc4y4         ones(size(dispbc4y4,1),1)  zeros(size(dispbc4y4,1),1)];
BC_Disp4     =    [BC_Dispx4;  BC_Dispy4];
dispbc4      =    [dispbc4x4;   dispbc4y4];
figure
plot(coord4(dispbc4,1),coord4(dispbc4,2),'b*' );
axis([0,6,0,6]);
title('bc-disp4')
%find point4 for force and pressure
point4number4   =    size(coord4,1);
point4         =    zeros(point4number4,1);
for i=1:point4number4
    if(abs(sqrt(coord4(i,1)^2+coord4(i,2)^2)-1)<eps)
        point4(i)=i;
        i=i+1;
    end
end
point4;
fpoint4=find(point4>eps);
figure
plot(coord4(fpoint4,1),coord4(fpoint4,2),'r+')
title('f-point4')
%find edge4
edge4=[];
sizesdconn4=length(sdConn4)
for i=1:sizesdconn4
    i;
    sdConn4{i};
    sizesdconn43=size(sdConn4{i},1);
    for  j=1:sizesdconn43
        lia=intersect(fpoint4,sdConn4{i}(j,:),'stable');
        if size(lia,1)==2
            edge4=[edge4;sdConn4{i}(j,:)];
        end
    end
end
edge4;
%calculate angle4
angle4=[];
fpoint4;
coord4(fpoint4,:);
x           =coord4(edge4,1);
y           =coord4(edge4,2);
angle4      =atan(y./x);
edge41      =edge4(:,1);
edge42      =edge4(:,2);
x1          =coord4(edge41,1);
y1          =coord4(edge41,2);
angle41     =atan(y1./x1);
x2          =coord4(edge42,1);
y2          =coord4(edge42,2);
angle42     =atan(y2./x2);
trac4       =[cos(angle41)   sin(angle41)   cos(angle42)   sin(angle42)]';
%add force or pressure
ndn4 = 2;
NDof4 = ndn4*size(coord4,1);
F4 = zeros(NDof4,1);
F4 = addSurfTraction(coord4, edge4, trac4, F4);
% solve
%ltx {\bf Solution of S-elements and assemblage of global stiffness}
%ltx {\bf and mass matrices}
[sdSln4, K4, M4] = SBFEMAssembly(coord4, sdConn4, sdsc4, mat);
%ltx {\bf Static solution of nodal displacements and forces}
[d4, F4] = SolverStatics(K4, BC_Disp4, F4);
d4;
d4(121);
d4(122);
d4(141);
d4(142);
%ltx {\bf Plot deformed mesh}
figure
opt = struct('MagnFct', 0.1, 'Undeformed','--k');
PlotDeformedMesh(d4, coord4, sdConn4, opt)
disp('Nodal displacements')
for ii = 1:length(coord4)
    fprintf('%5d %25.15e %25.15d\n',ii, d4(2*ii-1:2*ii))
end
sdSln4;
sdStrnMode4         =  SElementStrainMode2NodeEle( sdSln4 );%ltx integration constants
sdIntgConst4        =  SElementIntgConst( d4, sdSln4 );
%ltx displacements and strains at specified raidal coord4inate
d4;
sdsc4;
size(sdsc4,1);
nsd4=size(sdsc4,1);
for isd=1:nsd4
    xi = 1;    %ltx radial coord4inate
    [nodexy4{isd},dsp4{isd},strnNode4{isd},GPxy4{isd},strnEle4{isd}]=SElementInDispStrain(xi,sdSln4{isd},sdStrnMode4{isd},sdIntgConst4{isd});
end
for isd=1:nsd4
    ximid = 0.5;    %ltx radial coord4inate
    [nodexy4mid{isd},dsp4mid{isd},strnNode4mid{isd},GPxy4mid{isd},strnEle4mid{isd}]=SElementInDispStrain(ximid,sdSln4{isd},sdStrnMode4{isd},sdIntgConst4{isd});
end
for isd=1:nsd4
    xisdsc = 1e-15;    %ltx radial coord4inate
    [nodexy4sdsc{isd},dsp4sdsc{isd},strnNode4sdsc{isd},GPxy4sdsc{isd},strnEle4sdsc{isd}]=SElementInDispStrain(xisdsc,sdSln4{isd},sdStrnMode4{isd},sdIntgConst4{isd});
end
sdscstrain4=[];
for i=1:nsd4
    sdscstrain4(i,:)=strnNode4sdsc{i}(:,1)';
end
sdscstrain4;
%%
% % C = {1, 2, 3, 4, 5};  % 示例 cell 数组
% target = 61        % 要查找的数字
%
% % 示例 cell 数组，每个元素是一个矩阵
% C = sdConn4
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
%%
%   strnNode4{1707}
% Astrain=strnNode4{1707}(:,1)
% % Astress=mat.D * Astrain
% %  strnNode4{1118}
% Bstrain=strnNode4{1118}(:,3)
%
% Bstress=mat.D * Bstrain
magnfct44           = 0.1;
Umax4               = max(abs(d4));%ltx maximum displaceM4ent
Lmax4               = max(max(coord4)-min(coord4));
%ltx maximum dimension of domain
fct4                = magnfct44*Lmax4/Umax4;
%ltx factor to magnify the displaceM4ent ltx augment nodal coord4inates
dispmentforplot4    = fct4*(reshape(d4,2,[]))';
% stress recovery
[recoverystress4]  =spr(sdSln4,strnEle4,strnNode4,mat,nsd4,coord4,d4,GPxy4,strnEle4mid);
[sdscstressre4]=sprsdscstress(recoverystress4,sdSln4,coord4,nsd4,mat,strnNode4sdsc);
recoverystress4;
% error estimation use recovery stress
elementnumber4=nsd4;
[eta4,matupoly4,absepoly4,sdSln4]=errorestimation(sdIntgConst4,sdStrnMode4,strnNode4,strnEle4,elementnumber4,recoverystress4,mat,sdSln4,sdscstressre4,strnEle4mid,sdscstrain4);
% refine stragety 鎵惧埌闇€瑕佺粏鍖栫殑polygon锛岀劧鍚庡垽鏂槸鍑犺竟褰紝鐒跺悗缁嗗寲锛屽湪璁＄畻锛屽湪鍒ゆ柇锛岀洿鍒拌揪鍒扮粏鍖栬姹�
uall4=sum(matupoly4);
eall4=sum(absepoly4);
% reen4=sqrt(eall4)/sqrt(uall4)
%reen4=sqrt(0.5*eall4)/sqrt(0.5*uall4)
reen4=eall4/uall4
etamax4=max(absepoly4)
% eta4overline4                =sqrt(eall4^2/(uall4^2+eall4^2) );
% eta4max=0.1
% if eta4overline4 > eta4max%鍏ㄥ眬鍒ゆ柇
%     eM4            =eta4max*sqrt((uall4^2+eall4^2)/nsd4)     ;     %鍗曞厓骞冲潎鍏佽璇樊
%     ksi4           =cell(size(strnNode4,2),1);
%     for i=1:nsd4
%         ksi4{i}    = absepoly4(i)/ eM4;
%     end
%     marK4=[];
%     for i=1:nsd4
%         if ksi4{i} > 1
%             marK4(i)=i;
%         end
%     end
marK4=[]
for i=1:nsd4
    if absepoly4(i) > etamax4*gammy
        marK4(i)=i;
    end
end
marK4
size(marK4);
nonzeromarK4=nonzeros(marK4);
maxnodenumber4={};
for i=nonzeromarK4'
    i;
    disp(i);
    maxnodenumber4{i}=max(sdSln4{i}.node);
end

maxnodenumber4;
maxmaxnodenumber4= max(cell2mat(maxnodenumber4));
eleM4entsize4=[];
for i=nonzeromarK4'
    i;
    eleM4entsize4(i)=length(sdSln4{i}.node);
end
eleM4entsize4;
% eleM4entsize4(eleM4entsize4==0)=[];
length(eleM4entsize4);
addnodenumberstar4=[];
addnodenumberstar4(1)=maxmaxnodenumber4+1;
for j=2:length(nonzeromarK4)
    j;
    addnodenumberstar4(j)=addnodenumberstar4(j-1)+2*eleM4entsize4(nonzeromarK4(j-1));
end
addnodenumberstar4;
marK4;
% x=[2 8 4 9 3 ]
% y=[0 0 1 0 0 0 0 1 1 0 0 1 0 0 1]
% y(logical(y))=x
addnodenumberstar4y(logical(marK4))=addnodenumberstar4;
addnewnodecoord4={};
nonzeromarK4
for i=nonzeromarK4'
    [addnewnodecoord4{i}]=addnodecoord(sdSln4,nsd4,marK4,mat,i,coord4,sdConn4,sdsc4,addnodenumberstar4y);
end
addnewnodecoord4;
addnewnodecoord4M4=[];
for i=nonzeromarK4'
    i;
    AteM4p=[];
    AteM4p=addnewnodecoord4{i};
    addnewnodecoord4M4=[addnewnodecoord4M4;AteM4p];
end
addnewnodecoord4M4;
coord4;
teM4pdis=[];
for i=1:length(addnewnodecoord4M4)
    i;
    icoordx=addnewnodecoord4M4(i,2);
    icoordy=addnewnodecoord4M4(i,3);
    for j=i+1:length(addnewnodecoord4M4)
        j;
        jcoordx=addnewnodecoord4M4(j,2);
        jcoordy=addnewnodecoord4M4(j,3);
        teM4pdis=sqrt((jcoordx-icoordx)^2+(jcoordy-icoordy)^2);
        if teM4pdis<1e-5
            addnewnodecoord4M4(j,1)=addnewnodecoord4M4(i,1);
        end
    end
end
addnewnodecoord4M4;
addnewnodecoord4M4sort4=sortrows(addnewnodecoord4M4,1);
teM4pd=[];
for i=1:length(addnewnodecoord4M4sort4)
    i;
    inode=addnewnodecoord4M4sort4(i,1);
    for j=i+1:length(addnewnodecoord4M4sort4)
        j;
        jnode=addnewnodecoord4M4sort4(j,1);
        teM4pd=jnode-inode;
        if teM4pd<0.5
            addnewnodecoord4M4sort4(j,:)=[0 0 0];
        end
    end
end
addnewnodecoord4M4sort4;
size(addnewnodecoord4M4sort4);
addnewnodecoord4M4sort4(all(addnewnodecoord4M4sort4==0,2),:)=[];
size(addnewnodecoord4M4sort4);
addnewnodecoord4M4sort4;
addnewnodecoord4M4sort4u=addnewnodecoord4M4sort4;
coord4;
length(coord4);
addnewnodecoord4M4sort4urenode4=[addnewnodecoord4M4sort4u  (length(coord4)+1:length(coord4)+length(addnewnodecoord4M4sort4u))'];
addnewnodecoord4M4sort4urenode4;
addcoord4=[addnewnodecoord4M4sort4urenode4(:,end-2:end)];
addnewnodecoord4;
addnewnodepoly4={};
addnewnodecoord4refine4=[];
for i=nonzeromarK4'
    i;
    addnewnodecoord4refine4=addnewnodecoord4{i};
    for j=1:length(addnewnodecoord4refine4)
        j;
        jcoordx=addnewnodecoord4refine4(j,2);
        jcoordy=addnewnodecoord4refine4(j,3);
        tmepdisten=[];
        for k=1:length(addcoord4)
            kcoordx=addcoord4(k,1);
            kcoordy=addcoord4(k,2);
            tmepdisten=sqrt((kcoordx-jcoordx)^2+(kcoordy-jcoordy)^2);
            if tmepdisten<1e-5
                addnewnodecoord4refine4(j,1)=addcoord4(k,3);
                addnewnodepoly4{i}=addnewnodecoord4refine4;

            end
        end
    end
end
addnewnodepoly4;
sdConn4;
coord4;
sdsc4;
for i=nonzeromarK4'
    i;
    [coord4,sdConn4,sdsc4]=addnewpoly(addnewnodepoly4,sdSln4,nsd4,marK4,mat,i,maxmaxnodenumber4,coord4,sdConn4,sdsc4);
end
coord4;
sdConn4
sdsc4;
coord4=[coord4;addcoord4(:,1:2)];
%%鎵惧埌娌℃湁缁嗗寲鐨勫崟鍏冿紝缁欎粬缁嗗寲浜嗙殑杈逛笂澧炲姞涓€涓妭鐐�
%%娉ㄦ剰 鏈夌殑鍗曞厓浼氭湁澶氫釜杈硅缁嗗寲锛岄渶瑕佸鍔犲嚑涓偣
for i=1:length(sdConn4)
    i;
    ismember(nonzeromarK4,i);
    if ismember(nonzeromarK4,i)==0
        teM4peleM4ent =sdConn4{i};
        sdnode        =teM4peleM4ent(:,1);
        teM4pcoord4   =coord4(sdnode,:);
        for j= length(teM4peleM4ent):-1:1;
            j;
            teM4pnode       = teM4peleM4ent(j,:);
            teM4pnodecoord4 = [coord4(teM4pnode(1),1)+coord4(teM4pnode(2),1)   coord4(teM4pnode(1),2)+coord4(teM4pnode(2),2) ] /2 ;
            for k= 1:length(coord4)
                addnode        =k;
                addnodecoord4   =coord4(k,:);
                teM4pdisplment  =sqrt((coord4(k,1)-teM4pnodecoord4(1))^2+(coord4(k,2)-teM4pnodecoord4(2))^2);

                if  teM4pdisplment<1e-5
                    teM4peleM4ent(length(teM4peleM4ent)+1,:)=[0 0];
                    teM4peleM4ent(j+2:end,:)=teM4peleM4ent(j+1:end-1,:);
                    teM4peleM4ent(j,:)      =[teM4pnode(1) k];
                    teM4peleM4ent(j+1,:)    =[k           teM4pnode(2)];
                end
            end
        end
        teM4peleM4ent;
        sdConn4{i}=teM4peleM4ent;
    end
end
% elseif eta4overline4<=eta4max
%     disp('convergence')
% end
sdConn5=sdConn4;
sdsc5=sdsc4;
coord5=coord4;
% figure
% opt        = struct('LineSpec','-k', 'sdSC',sdsc5);
% for i=1:size(sdConn5,1)
%     sdsc5(i,:)
%     sdsc5(i,1)
%     sdsc5(i,2)
%     str = i
%     text(sdsc5(i,1),sdsc5(i,2),num2str(str))
% end
% PlotSBFEMesh(coord5, sdConn5, opt);
% title('MESH5');
figure
opt        = struct('LineSpec','-k', 'sdSC',sdsc5);
PlotSBFEMesh(coord5, sdConn5, opt);
title('MESH5');
% figure
% plot(coord5(:,1),coord5(:,2),'b+')
% axis equal
% title('point5')



