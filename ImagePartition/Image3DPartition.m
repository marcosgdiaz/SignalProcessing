function [ ] = Image3DPartition( x, y, z, rho, umbral, volumen, volumen2 )
binaryrho = rho > umbral;
filledrho = imfill(binaryrho,'holes');
CC=bwconncomp(filledrho,26);
prop = regionprops(CC, 'FilledArea','PixelIdxList');
allAreas = [prop.FilledArea]
allowableAreaIndexes = allAreas > volumen & allAreas < volumen2; 
keeperIndexes = find(allowableAreaIndexes);

for i=keeperIndexes
    img=zeros(size(rho));
    aux=prop(i).PixelIdxList;
    img(aux)=true;
    cdata = rand()*ones(size(img));
    p=patch(isosurface(x,y,z,img,0));
    isonormals(x,y,z,img,p)
    isocolors(x,y,z,cdata,p)
    p.FaceColor = 'interp';
    p.EdgeColor = 'none';
    grid on
    hold all
end
%keeperBlobsImage = ismember(filledrho, keeperIndexes);
view(150,30)
set(gca,'Zdir','reverse')
xlabel('x(m)');
ylabel('y(m)');
zlabel('z(m)')
title('3-D binary image')
end

