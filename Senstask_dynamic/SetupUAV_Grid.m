function Grid=SetupUAV_Grid(Xlim,dx)


[meshX,meshY]=meshgrid(Xlim(1):dx:Xlim(2));
%             XY=[Y(:),X(:)];
XY=zeros(size(meshX,1)*size(meshX,2),2);
XYindex=zeros(size(meshX,1)*size(meshX,2),1);
XYgridpos=zeros(size(meshX,1)*size(meshX,2),2);
k=1;
for i=1:1:size(meshX,1)
    for j=1:1:size(meshX,2)
        XY(k,:)=[meshY(i,j),meshX(i,j)];
        XYindex(k)=k;
        XYgridpos(k,:)=[i,j];
        k=k+1;
    end
end
%             XYindex=[1:1:size(XY,1)]';
% [XY,XYindex]
ng=length(XYindex);
XYADJindex=cell(ng,1);

for i=1:1:length(XYindex)
    xy=XY(i,:);
    p=[];
    
    r=xy+[dx,0];
    if r(1)>Xlim(2) || r(2) >Xlim(2) || r(1)<Xlim(1) || r(2)<Xlim(1)
        
    else
        [~,ind]=min(sqrt(sum((repmat(r,ng,1)-XY).^2,2)));
        p=horzcat(p,ind(1));
    end
    
    r=xy+[0,dx];
    if r(1)>Xlim(2) || r(2) >Xlim(2) || r(1)<Xlim(1) || r(2)<Xlim(1)
        
    else
        [~,ind]=min(sqrt(sum((repmat(r,ng,1)-XY).^2,2)));
        p=horzcat(p,ind(1));
    end
    
    r=xy+[-dx,0];
    if r(1)>Xlim(2) || r(2) >Xlim(2) || r(1)<Xlim(1) || r(2)<Xlim(1)
        
    else
        [~,ind]=min(sqrt(sum((repmat(r,ng,1)-XY).^2,2)));
        p=horzcat(p,ind(1));
    end
    
    r=xy+[0,-dx];
    if r(1)>Xlim(2) || r(2) >Xlim(2) || r(1)<Xlim(1) || r(2)<Xlim(1)
        
    else
        [~,ind]=min(sqrt(sum((repmat(r,ng,1)-XY).^2,2)));
        p=horzcat(p,ind(1));
    end
    
    XYADJindex{i}=sort(p);
end

Grid.meshX=meshX;
Grid.meshY=meshY;
Grid.XY=XY;
Grid.XYindex=XYindex;
Grid.XYADJindex=XYADJindex;
Grid.XYgridpos=XYgridpos;
Grid.Xlim = Xlim;
Grid.dx = dx;


end


