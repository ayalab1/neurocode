function [av,st,pts]=brainPad(av,st,pts)
%% EDIT THIS LINE TO THE LOCATION OF DOWNLOADED ANNOTATED VOLUME AND
% STRUCTURE TREE, OR ELSE PASS IN THE FOLDER PATH AS INPUT
%volumePath = '/Users/timothyspellman/Dropbox/MATLAB/brainPad';
%
% INPUTS
% av:      10um annotation volume, [AP,DV,ML], Allen CCF, .npy file -- 
%          pass this input for faster loading, or load from within script
% st:      structure tree, .csv file
% pts:  AP/ML/DV annotated coordinates, bregma-centered for AP/ML,
%          brain-surface-zeroed for DV
%
% Requirements
% npy-matlab: https://github.com/kwikteam/npy-matlab
% (Copyright (c) 2015, npy-matlab developers, All rights reserved)
% annotation_volume_10um_by_index.npy, structure_tree_safe_2017.csv: http://data.cortexlab.net/allenCCF
% 
% References
% https://www.biorxiv.org/content/early/2018/10/19/447995
% This code incorporates elements of the allenCCF package from Cortex Lab and Kwik Team:
% https://github.com/cortex-lab/allenCCF
% http://help.brain-map.org/download/attachments/2818169/MouseCCF.pdf
%
%% Tim Spellman, 2021

viewAxis='coronal';
x=534; %%Bregma in CCF space is 534,570
t={};
mouse_pos=[400,570];

if nargin<2
    av=readNPY([volumePath filesep 'annotation_volume_10um_by_index.npy']);
    st=openStructure([volumePath filesep 'structure_tree_safe_2017.csv']);
else
    if ischar(av);
        volumePath=av;
    end
end
assignin('base','av',av);
assignin('base','st',st);

fprintf('\n\n');
fprintf('   [left / right arrows] = move forward/back \n')
fprintf('   [control] = jump forward/back by 100um \n')
fprintf('   [c/s/h] = change view axis \n')
fprintf('   [click] = add/remove points \n')
fprintf('   [d] = connect-the-dots on/off \n\n')


depthChart=zeros(size(av,1),size(av,3));
for xx=1:size(av,2)/2
   depthChart(find(squeeze(av(:,xx,:))~=1  & depthChart==0 ))=xx;
end;

%%initialize coordinate field or convert from bregma xyz(um) to voxels
if nargin<3
    pts=zeros(3,0);
else
    pts(1,:,:)=round((pts(1,:,:)-5340)./-10);
    pts(3,:,:)=round((pts(3,:,:)+5700)./10);
    pts(2,:,:)=round(pts(2,:,:)./-10+depthChart(sub2ind(size(depthChart),pts(1,:,:),pts(3,:,:))));
end
d=0;
f=figure('keyPressFcn',@keyPress,'WindowButtonMotionFcn',@mouseOver,...
    'WindowButtonDownFcn',@leftClick);
set(f,'Pointer','crosshair');
guidata(f,[{x},{av},{viewAxis},{st},{t},{depthChart},{mouse_pos},{pts},{d}]);
refreshFig(f)
end

function refreshFig(f)
gdata=guidata(f);
x=gdata{1};
av=gdata{2};
viewAxis=gdata{3};
st=gdata{4};
t=gdata{5};
depthChart=gdata{6};
mouse_pos=gdata{7};
pts=gdata{8};
d=gdata{9};
if contains(viewAxis,'coronal')
    imav=squeeze(av(x,:,:));
    imav=~(round(conv2(imav,ones(3)./9,'same'))~=imav);
    imshow(imav+0.7);
    ax1=gca;colormap(ax1,gray);
    hold on;
    crit=find((pts(1,:)>=x-3 & pts(1,:)<=x+3) & pts(1,:)~=x);
    sc(1)=scatter(pts(3,crit),pts(2,crit),'ok','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0,'MarkerFaceColor','r'); 
    crit=find((pts(1,:)==x-2 | pts(1,:)==x+2) & pts(1,:)~=x);
    sc(2)=scatter(pts(3,crit),pts(2,crit),'ok','MarkerFaceAlpha',0.4,'MarkerEdgeAlpha',0,'MarkerFaceColor','r');
    crit=find((pts(1,:)==x-1 | pts(1,:)==x+1) & pts(1,:)~=x);
    sc(3)=scatter(pts(3,crit),pts(2,crit),'ok','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0,'MarkerFaceColor','r');
    crit=find((pts(1,:)==x));
    sc(4)=scatter(pts(3,crit),pts(2,crit),'ok','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0,'MarkerFaceColor','r');
    if d 
        for q=1:4;sc(q).SizeData=5;end
      crit=find(pts(1,:)>=x-3 & pts(1,:)<=x+3);
      try;ch=convhull(pts(3,crit),pts(2,crit));catch;ch=[];end
      gcf;hold on;ptch=patch(pts(3,crit(ch)),pts(2,crit(ch)),'r','EdgeAlpha',0,'FaceAlpha',0.2);  
    end
    hold off;
    
elseif contains(viewAxis,'sagittal')
    imav=squeeze(av(:,:,x))';
    imav=~(round(conv2(imav,ones(3)./9,'same'))~=imav);
    imshow(imav+0.7);
    ax1=gca;colormap(ax1,gray);
    hold on;
    crit=find((pts(3,:)==x-3 | pts(3,:)==x+3) & pts(3,:)~=x);
    sc(1)=scatter(pts(1,crit),pts(2,crit),'ok','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0,'MarkerFaceColor','r');
    crit=find((pts(3,:)==x-2 | pts(3,:)==x+2) & pts(3,:)~=x);
    sc(2)=scatter(pts(1,crit),pts(2,crit),'ok','MarkerFaceAlpha',0.4,'MarkerEdgeAlpha',0,'MarkerFaceColor','r');
    crit=find((pts(3,:)==x-1 | pts(3,:)==x+1) & pts(3,:)~=x);
    sc(3)=scatter(pts(1,crit),pts(2,crit),'ok','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0,'MarkerFaceColor','r');
    crit=find((pts(3,:)==x));
    sc(4)=scatter(pts(1,crit),pts(2,crit),'ok','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0,'MarkerFaceColor','r');
    if d 
        for q=1:4;sc(q).SizeData=5;end
        crit=find(pts(3,:)>=x-3 & pts(3,:)<=x+3);
        try;ch=convhull(pts(1,crit),pts(2,crit));catch;ch=[];end
        gcf;hold on;ptch=patch(pts(1,crit(ch)),pts(2,crit(ch)),'r','EdgeAlpha',0,'FaceAlpha',0.2);  
    end
    hold off;
elseif contains(viewAxis,'horizontal')
    imav=squeeze(av(:,x,:))';
    imav=~(round(conv2(imav,ones(3)./9,'same'))~=imav);
    imshow(imav+0.7);
    ax1=gca;colormap(ax1,gray);
    hold on;
    crit=find((pts(2,:)==x-3 | pts(2,:)==x+3) & pts(2,:)~=x);
    sc(1)=scatter(pts(1,crit),pts(3,crit),'ok','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0,'MarkerFaceColor','r');
    crit=find((pts(2,:)==x-2 | pts(2,:)==x+2) & pts(2,:)~=x);
    sc(2)=scatter(pts(1,crit),pts(3,crit),'ok','MarkerFaceAlpha',0.4,'MarkerEdgeAlpha',0,'MarkerFaceColor','r');
    crit=find((pts(2,:)==x-1 | pts(2,:)==x+1) & pts(2,:)~=x);
    sc(3)=scatter(pts(1,crit),pts(3,crit),'ok','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0,'MarkerFaceColor','r');
    crit=find((pts(2,:)==x));
    sc(4)=scatter(pts(1,crit),pts(3,crit),'ok','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0,'MarkerFaceColor','r');
    if d 
        for q=1:4;sc(q).SizeData=5;end
        crit=find(pts(2,:)>=x-3 & pts(2,:)<=x+3);
        try;ch=convhull(pts(1,crit),pts(3,crit));catch;ch=[];end
        gcf;hold on;ptch=patch(pts(1,crit(ch)),pts(3,crit(ch)),'r','EdgeAlpha',0,'FaceAlpha',0.2);  
    end
    hold off;
end
gcf;
if contains(viewAxis,'coronal')
    text(25,50,['AP: ',num2str(x*-10+5340),'μm'])
    if ~isempty(t); 
        text(25,75,t{2});
        text(25,100,t{3});
    else
    text(25,75,['ML: ',num2str(mouse_pos(2)*10-5700),'μm']);
    text(25,100,['DV: ',num2str(mouse_pos(1)*10+depthChart(mouse_pos(2))*10),'μm']);
    end
elseif contains(viewAxis,'sagittal')
    text(25,85,['ML: ',num2str(x*10-5700),'μm'])
    if ~isempty(t);
        text(25,55,t{1});
        text(25,115,t{3});
    else
    text(25,55,['AP: ',num2str(mouse_pos(2)*10+5340),'μm']);
    text(25,115,['DV: ',num2str(mouse_pos(1)*10+depthChart(mouse_pos(2))*10),'μm']);
    end
elseif contains(viewAxis,'horizontal')
    text(25,115,['DV: ',num2str(x*-10+340),'μm'])
    if ~isempty(t);
        text(25,55,t{1});
        text(25,85,t{2});
    else
    text(25,55,['AP: ',num2str(mouse_pos(2)*10+5340),'μm']);
    text(25,85,['ML: ',num2str(mouse_pos(1)*10-5700),'μm']);
    end
end
    
hold off
end

function keyPress(f,event)
gdata=guidata(f);
x=gdata{1};
av=gdata{2};
viewAxis=gdata{3};
st=gdata{4};
t=gdata{5};
depthChart=gdata{6};
mouse_pos=gdata{7};
pts=gdata{8};
d=gdata{9};
switch event.Key
       case 'leftarrow'
        if contains(event.Modifier,'control')
        try x=x-10;catch display('First slice!');end
        else try x=x-1;catch display('First slice!');end
        end
        guidata(f,[{x},{av},{viewAxis},{st},{t},{depthChart},{mouse_pos},{pts},{d}]);
        refreshFig(f);
       case 'rightarrow'
        if contains(event.Modifier,'control')
        try x=x+10;catch display('Last slice!');end
        else try x=x+1;catch display('First slice!');end
        end
        guidata(f,[{x},{av},{viewAxis},{st},{t},{depthChart},{mouse_pos},{pts},{d}]);
        refreshFig(f);
       case 'c'
           if contains(viewAxis,'sagittal')
               x=mouse_pos(1);
           elseif contains(viewAxis,'horizontal')
               x=mouse_pos(1);
           end
           viewAxis='coronal';
           guidata(f,[{x},{av},{viewAxis},{st},{t},{depthChart},{mouse_pos},{pts},{d}]);
            refreshFig(f);
        case 's'
            if contains(viewAxis,'coronal')
               x=mouse_pos(1);
           elseif contains(viewAxis,'horizontal')
               x=mouse_pos(2);
           end
           viewAxis='sagittal';
           guidata(f,[{x},{av},{viewAxis},{st},{t},{depthChart},{mouse_pos},{pts},{d}]);
           refreshFig(f);
        case 'h'
            if contains(viewAxis,'coronal')
               x=mouse_pos(2);
           elseif contains(viewAxis,'sagittal')
               x=mouse_pos(2);
           end
           viewAxis='horizontal';
           guidata(f,[{x},{av},{viewAxis},{st},{t},{depthChart},{mouse_pos},{pts},{d}]);
           refreshFig(f);
        case 'd'
            if d==0;d=1;else d=0;end
            guidata(f,[{x},{av},{viewAxis},{st},{t},{depthChart},{mouse_pos},{pts},{d}]);
            refreshFig(f);
            
end;end

function mouseOver(f,eventdata)
gdata=guidata(f);
x=gdata{1};
av=gdata{2};
viewAxis=gdata{3};
st=gdata{4};
t=gdata{5};
depthChart=gdata{6};
mouse_pos=gdata{7};
pts=gdata{8};
d=gdata{9};

ax=findall(f,'Type','Axes');
if contains(viewAxis,'coronal')
im=squeeze(av(x,:,:));
dc=squeeze(depthChart(x,:));
elseif contains(viewAxis,'sagittal')
im=squeeze(av(:,:,x))';
dc=squeeze(depthChart(:,x))';
elseif contains(viewAxis,'horizontal')
im=squeeze(av(:,x,:))';
dc=depthChart';
end
mouse_position=get(ax(1),'CurrentPoint');
mX=round(mouse_position(1,1));
mY=round(mouse_position(1,2));

if contains(viewAxis,'coronal') &...
        mX>0 & mX<=length(dc)
t(1)={['AP: ',num2str(x*-10+5340),'µm']};
t(2)={['ML: ',num2str(mX*10-5700),'µm']};
t(3)={['DV: ',num2str(mY*-10+dc(mX)*10),'µm']};
elseif contains(viewAxis,'sagittal') &...
        mX>0 & mX<=length(dc)
t(1)={['AP: ',num2str(mX*-10+5340),'µm']};
t(2)={['ML: ',num2str(x*10-5700),'µm']};
t(3)={['DV: ',num2str(mY*-10+dc(mX)*10),'µm']};
elseif contains(viewAxis,'horizontal') &...
        mX>0 & mX<=size(dc,2) & mY>0 & mY<=size(dc,1);
t(1)={['AP: ',num2str(mX*-10+5340),'µm']};
t(2)={['ML: ',num2str(mY*10-5700),'µm']};
t(3)={['DV: ',num2str(x*-10+dc(mY,mX)*10),'µm']};
end

if ~ismember(mX,1:size(im,2)) |...
        ~ismember(mY,1:size(im,1))
    return
end
currAV=im(mY,mX);
if currAV==1
    return
end
mouse_pos=[mX,mY];

areaName=st.safe_name(currAV);
gcf;
fa=findall(f,'Type','Text');
if ~isempty(t);
    delete(fa)
    if contains(viewAxis,'coronal')
    text(25,50,t(1));
    text(25,75,t(2));
    text(25,100,t(3));
    else
    text(25,55,t(1));
    text(25,85,t(2));
    text(25,115,t(3));
    end
end
text(25,25,areaName);
guidata(f,[{x},{av},{viewAxis},{st},{t},{depthChart},{mouse_pos},{pts},{d}]);
end

function leftClick(f,eventdata)
gdata=guidata(f);
x=gdata{1};
av=gdata{2};
viewAxis=gdata{3};
st=gdata{4};
t=gdata{5};
depthChart=gdata{6};
mouse_pos=gdata{7};
pts=gdata{8};
d=gdata{9};

cp=get(gca,'CurrentPoint'); cp=round(cp(1,1:2));
if contains(viewAxis,'coronal')
    ptsNew=[x;cp(2);cp(1)];
elseif contains(viewAxis,'sagittal')
    ptsNew=[cp(1);cp(2);x];
elseif contains(viewAxis,'horizontal')
    ptsNew=[cp(1);x;cp(2)];
end
 
rm=[];
for q=1:size(pts,2)
    if nthroot( (ptsNew(1)-pts(1,q))^2 + (ptsNew(2)-pts(2,q))^2 + (ptsNew(3)-pts(3,q))^2, 3)<2;
        rm=[rm,q];
end;end
if isempty(rm);pts=[pts,ptsNew];else pts(:,rm)=[];end
          
ptsOut=zeros(size(pts));
    ptsOut(1,:,:)=pts(1,:,:).*-10+5340;
    ptsOut(3,:,:)=pts(3,:,:).*10-5700;
    ptsOut(2,:,:)=pts(2,:,:).*-10+depthChart(sub2ind(size(depthChart),pts(1,:,:),pts(3,:,:))).*10;
assignin('base','pts',ptsOut);
guidata(f,[{x},{av},{viewAxis},{st},{t},{depthChart},{mouse_pos},{pts},{d}]);
refreshFig(f);
end

function structure=openStructure(fn)
[~, fnBase] = fileparts(fn);
fid = fopen(fn, 'r');
    titles = textscan(fid,repmat('%s',1,21),1,'delimiter', ',');
    titles = cellfun(@(x)x{1}, titles, 'uni', false);  
    data = textscan(fid, ['%d%d%s%s'... % 'id'    'atlas_id'    'name'    'acronym'
                          '%s%d%d%d'... % 'st_level'    'ontology_id'    'hemisphere_id'    'weight'
                          '%d%d%d%d'... % 'parent_structure_id'    'depth'    'graph_id'     'graph_order'
                          '%s%s%d%s'... % 'structure_id_path'    'color_hex_triplet' neuro_name_structure_id neuro_name_structure_id_path
                          '%s%d%d%d'... % 'failed'    'sphinx_id' structure_name_facet failed_facet
                          '%s'],'delimiter',','); % safe_name
    titles = ['index' titles];
    data = [[0:numel(data{1})-1]' data];    
structure=table(data{:},'VariableNames',titles);
fclose(fid);
end

