% Generate points that fill space evenly but randomly, using 2d poisson disc
% Rhodri Cusack, Brain and Mind Institute, Western University, Canada, July 2013
% www.cusacklab.org  cusacklabcomputing.blogspot.ca rhodri@cusacklab.org
% 
% Algorithm with thanks from
% http://devmag.org.za/2009/05/03/poisson-disk-sampling/
%
% sz=[width, height] of space to be filled
% min_dist=min separation of points
% newpointscount= higher number gives higher quality (fewer gaps)
%
% Example useage:
% spoints=generate_poisson_2d([100 100 ],10,20);

function [samplepoints]=generate_poisson_2d(sz,min_dist,newpointscount)

cellsize=min_dist/sqrt(2);
grid=cell(ceil(sz(1)/cellsize),ceil(sz(2)/cellsize));
proclist=[];
samplepoints=[];

% Random start
firstpoint=ceil(sz.*rand(1,2));

% This will be a queue with points pulled randomly from it
proclist=[proclist; firstpoint];

% Output...
samplepoints=[samplepoints; firstpoint];

% Grid - see algorithm from devmag above
gridpoint=imageToGrid(firstpoint,cellsize);
grid{gridpoint(1),gridpoint(2)}=firstpoint;

while ~isempty(proclist)
    randrow=ceil(rand(1)*size(proclist,1));
    point=proclist(randrow,:);
    proclist(randrow,:)=[];
    
    for i=1:newpointscount
        newpoint=generateRandomPointsAround(point, min_dist);
        if inRectangle(newpoint,sz) && ~inNeighbourhood(grid, newpoint, min_dist,cellsize)
            proclist=[proclist; newpoint];
            samplepoints=[samplepoints; newpoint];
            gridpoint=imageToGrid(newpoint,cellsize);
            grid{gridpoint(1),gridpoint(2)}=newpoint;
        end;
    end;
end;


figure(10);
scatter(samplepoints(:,1),samplepoints(:,2));


end

function [gpoint]=imageToGrid(point,cellsize)
gpoint=ceil(point/cellsize);
end

function [newpoint]=generateRandomPointsAround(point,min_dist)
[x y z]=sph2cart(2*pi*rand(1),0,min_dist*(rand(1)+1));
newpoint=point+[x y];
end

% Is there another point already nearby. 
function [isin]=inNeighbourhood(grid,point,min_dist,cellsize)
gridsz=size(grid);
% Where does this point belong in the grid
gridpoint=imageToGrid(point,cellsize);
% only check neighbours -2<delta<2 in each dim arount "gridpoint"
[ox oy]=meshgrid(-2:2,-2:2); 
c=repmat(gridpoint,[size(ox(:),1) 1])+[ox(:) oy(:)];
% Reject any putative neighbours that are out of bounds?
c(any(c<1,2) | c(:,1)>gridsz(1) | c(:,2)>gridsz(2),:)=[];
% Reject any putative neighbours without coordinates
c(isempty(cat(1,grid{sub2ind(gridsz,c(:,1),c(:,2))})),:)=[];
% Get points from grid neighbours 
neighbour_points=cat(1,grid{sub2ind(gridsz,c(:,1),c(:,2))});
% Any closeby?
if ~isempty(neighbour_points)
    dists=sqrt(sum((neighbour_points-repmat(point,[size(neighbour_points,1) 1])).^2,2));
    isin=any(dists<min_dist);
else
    isin=false;
end;
end

% Is point in rectangle specified by sz?
function [isin]=inRectangle(point,sz)
isin=all(point>1) && all(point<=sz);
end


