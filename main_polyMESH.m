% main_polyMESH.m
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: Routine for the generation of finite element meshes of 
%          polycrystalline microstructures with fully bonded grains
%
%
% This code is part of the following article. If used for other research,
% please acknowledge this contribution by citing the following reference: 
%
% M. Paggi, M. Corrado, J. Reinoso (2017) Fracture of solar-grade anisotropic
% polycrystalline Silicon: A combined phase field - cohesive zone model approach,
% Computer Methods in Applied Mechanics and Engineering, 
% Volume 330, 1 March 2018, Pages 123-148,
% https://doi.org/10.1016/j.cma.2017.10.021  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input data
dim_x=20;           % dimension of the specimen in the x-direction
dim_y=20;           % dimension of the specimen in the y-direction
nb_grains=10;       % number of grains to be generated in the specimen
quadrangle=true;    % true: the specimen will be meshed with quadrangular elements
                    % false: the specimen will be meshed with triangular elements
sf=1;               % scale factor, to rescale the specimen size
lc=2;               % approximate finite element size, for Gmsh


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of the polygrain microstructure with the Voronoi algorithm
rand('state',1);    % the second parameter is a kind of seed for the random
                    % set of points generated by rand()
x=rand(1,nb_grains)*dim_x;
y=rand(1,nb_grains)*dim_y;
[vx,vy] = voronoi(x,y,[]);
[v,c]=voronoin([x(:) y(:)]);
plot(x,y,'r+',vx,vy,'b-'); axis equal
xlim([0 dim_x]);
ylim([0 dim_y]);

% Construction of the connectivity matrix of the grains (connect_gr)
connect_gr=zeros(nb_grains,100);
for i=1:length(c)
    prov=c{i,1};
    vertices_per_grains(i)=length(prov);    % vertices_per_grains contains the number of vertices for each grain
    for j=1:length(prov)
        connect_gr(i,j)=prov(j);            % connectivity matrix of the grains
    end
end

nb_segments=length(vx);  % total number of segments
nb_nodes=length(v);      % total number of vertices

% The Voronoi algorithm generates grains that are extended beyond the edges
% of the specimen. They have to be modified by limiting their size in
% correspondence of the edges of the specimen.

% The segments of the grains that intersect the edges of the specimen are
% identified
k=0;
for i=1:nb_segments
    if ((vx(1,i)<0) | (vx(1,i)>dim_x) |(vx(2,i)<0) | (vx(2,i)>dim_x) | (vy(1,i)<0) | (vy(1,i)>dim_y) | (vy(2,i)<0) | (vy(2,i)>dim_y))
        if ((((vx(1,i)>0) & (vx(1,i)<dim_x)) | ((vx(2,i)>0) & (vx(2,i)<dim_x))) & (((vy(1,i)>0) & (vy(1,i)<dim_y)) | ((vy(2,i)>0) & (vy(2,i)<dim_y))))
            k=k+1;
            out(k)=i;   % segments intersecting the edges of the specimen
        end
    end
end

% Nodes 1 e 2 of the segments intersecting the edges of the specimen
node1=zeros(k,1);
node2=zeros(k,1);

% The nodes that are outside the domain of the specimen are identified
for i=1:k
    for j=1:nb_nodes
        if (abs(vx(1,out(i))-v(j,1))<0.00001 & abs(vy(1,out(i))-v(j,2))<0.00001)
            node1(i)=j;
        end
        if (abs(vx(2,out(i))-v(j,1))<0.00001 & abs(vy(2,out(i))-v(j,2))<0.00001)
            node2(i)=j;
        end
    end
    if node1(i)==0
        node1(i)=1;
    end
    if node2(i)==0
        node2(i)=1;
    end
end

% For each segment that intersects the domain the adiacent grains are
% identified
for i=1:k
    w=0;
    for j=1:nb_grains
        flag=0;
        for h=1:10
            if (flag==0 & node1(i)==connect_gr(j,h))
                for l=1:10
                    if (node2(i)==connect_gr(j,l))
                        w=w+1;
                        gr(w)=j;
                        flag=1;
                    end
                end
            end
        end
    end
    
    % For each intersecting segment the edge of intersection is determined.
    % The id_numbers of the edges of the specimen are the followings:
    %         ____3___
    %        |        |
    %      4 |        | 2
    %        |        |
    %        |________|
    %             1 
    edge=0;          
    A=[(vx(1,out(i))) (vx(2,out(i)))];
    if max(A)>dim_x
        xint=dim_x;
        yint=vy(1,out(i))+(xint-vx(1,out(i)))*(vy(2,out(i))-vy(1,out(i)))/(vx(2,out(i))-vx(1,out(i)));
        if (yint>0 & yint<dim_y)
            edge=2;
        end
    elseif min(A)<0
        xint=0;
        yint=vy(1,out(i))+(xint-vx(1,out(i)))*(vy(2,out(i))-vy(1,out(i)))/(vx(2,out(i))-vx(1,out(i)));
        if (yint>0 & yint<dim_y)
            edge=4;
        end
    end
    B=[(vy(1,out(i))) (vy(2,out(i)))];
    if max(B)>dim_y
        yint=dim_y;
        xint=vx(1,out(i))+(yint-vy(1,out(i)))*(vx(2,out(i))-vx(1,out(i)))/(vy(2,out(i))-vy(1,out(i)));
        if (xint>0 & xint<dim_x)
            edge=3;
        end
    elseif min(B)<0
        yint=0;
        xint=vx(1,out(i))+(yint-vy(1,out(i)))*(vx(2,out(i))-vx(1,out(i)))/(vy(2,out(i))-vy(1,out(i)));
        if (xint>0 & xint<dim_x)
            edge=1;
        end
    end
    
    % The coordinates of the intersection points are determined
    if edge==1
        yint=0;
        xint=vx(1,out(i))+(yint-vy(1,out(i)))*(vx(2,out(i))-vx(1,out(i)))/(vy(2,out(i))-vy(1,out(i)));
    elseif edge==2
        xint=dim_x;
        yint=vy(1,out(i))+(xint-vx(1,out(i)))*(vy(2,out(i))-vy(1,out(i)))/(vx(2,out(i))-vx(1,out(i)));
    elseif edge==3
        yint=dim_y;
        xint=vx(1,out(i))+(yint-vy(1,out(i)))*(vx(2,out(i))-vx(1,out(i)))/(vy(2,out(i))-vy(1,out(i)));
    elseif edge==4
        xint=0;
        yint=vy(1,out(i))+(xint-vx(1,out(i)))*(vy(2,out(i))-vy(1,out(i)))/(vx(2,out(i))-vx(1,out(i)));
    end
    
    % The connectivity matrix of the grains is modified by adding the points 
    % of intersection between the original segments and the edges of the
    % specimen
    for h=1:w
        connect_gr(gr(h),vertices_per_grains(gr(h))+1)=nb_nodes+1;
        vertices_per_grains(gr(h))=vertices_per_grains(gr(h))+1;
    end
    nb_nodes=nb_nodes+1;
    % the coordinates of the new vertex are inserted in the vector v
    v(nb_nodes,1)=xint;
    v(nb_nodes,2)=yint;
end

corners=[0 0; dim_x 0; dim_x dim_y; 0 dim_y];
for i=1:4
    dist_x=0;
    dist_y=0;
    min_dist_x=dim_x;
    min_dist_y=dim_y;
    for j=1:nb_nodes
        if (v(j,2)==corners(i,2))
            dist_x=abs(corners(i,1)-v(j,1));
            if (dist_x < min_dist_x)
                min_dist_x=dist_x;
                idx=j;
            end
        end
        
        if v(j,1)==corners(i,1)
            dist_y=abs(corners(i,2)-v(j,2));
            if dist_y < min_dist_y
                min_dist_y=dist_y;
                idy=j;
            end
        end
    end

    for j=1:nb_grains
        if (find(idx-connect_gr(j,1:vertices_per_grains(j))==0))
            if (find(idy-connect_gr(j,1:vertices_per_grains(j))==0))
                vertices_per_grains(j)=vertices_per_grains(j)+1;
                connect_gr(j,vertices_per_grains(j))=nb_nodes+1;
            end
        end
    end
nb_nodes=nb_nodes+1;
v(nb_nodes,:)=corners(i,:);
end

% The vertices outside the specimen edges are deleted from the
% connectivity matrix of the grains
for i=1:nb_grains
    k=1;
    for j=1:100
        if connect_gr(i,k)~=0
            if (v(connect_gr(i,k),1)<0 | v(connect_gr(i,k),1)>dim_x | v(connect_gr(i,k),2)<0 | v(connect_gr(i,k),2)>dim_y)
                for h=k:19
                    connect_gr(i,h)=connect_gr(i,h+1);
                end
                vertices_per_grains(i)=vertices_per_grains(i)-1;
                k=k;
            else
                k=k+1;
            end
        end
    end
end

% The vertices of the grains are ordered counterclockwise
for i=1:nb_grains
    centre(i,1)=x(i);
    centre(i,2)=y(i);
    for j=1:vertices_per_grains(i)
        deltay=(v(connect_gr(i,j),2)-centre(i,2));
        deltax=(v(connect_gr(i,j),1)-centre(i,1));
        tang(j)=deltay/deltax;
        if (deltay>0 & deltax>0)
            teta(j)=atan(tang(j));
        elseif (deltay>0 & deltax<0)|(deltay<0 & deltax<0)
            teta(j)=atan(tang(j))+pi;
        elseif (deltay<0 & deltax>0)
            teta(j)=atan(tang(j))+2*pi;
        end
    end
    prov=zeros(1,100);
    for j=1:vertices_per_grains(i)
        [C,index]=min(teta);
        prov(j)=connect_gr(i,index);
        teta(index)=10;
    end
    connect_gr(i,:)=prov(:);
    connect_gr(i,vertices_per_grains(i)+1)=connect_gr(i,1);
end

% The nodes outside the domain of the specimen are deleted from the vector
% of the nodal coordinates (coord); the remaining nodes are re-numbered and
% the connectivity matrix is modified accordingly
nb_nodes=length(v);
ind=1;
for i=1:length(v)
    nb_nodes=nb_nodes;
    if (i==1 | v(i,1)<0 | v(i,1)>dim_x | v(i,2)<0 | v(i,2)>dim_y)
        for j=ind:nb_nodes-1
            coord(j,:)=v(j+i-ind+1,:);
            new_id(j+i-ind+1,1)=j;
        end
        new_id(i)=0;
        nb_nodes=nb_nodes-1;
    else
        ind=ind+1;
    end
end

for i=1:nb_grains
    for j=1:vertices_per_grains(i)+1
        connect_gr(i,j)=new_id(connect_gr(i,j));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of the geometrical entities for Gmsh

% The model is re-scaled (if necessary)
coord=coord*sf;
nb_nodes=length(coord);

% Lines and lineloops are constructed
nline=0;
for i=1:nb_grains
    for j=1:vertices_per_grains(i)
        node1=connect_gr(i,j);
        node2=connect_gr(i,j+1);
		sign=0;
		for k=1:nline
		    if(mline(k,1)==node1)
		        if (mline(k,2)==node2)
		            sign=+1;
		            numline=k;
		        end
		    end
		    if(mline(k,1)==node2)
		        if (mline(k,2)==node1)
		            sign=-1;
		            numline=k;
		        end
		    end
		end
		if (sign==0)
		    nline=nline+1;
            mline(nline,1)=node1;
            mline(nline,2)=node2;
            numline=nline;
            sign=+1;
        end
        lineloop(i,j)=numline*sign;
    end
end

% The lines constituting the grain boundaries, except those coincident with 
% one of the specimen edges, are stored in a separate vector ('int_lines')
k=0;
h=0;
for i=1:length(mline)
    if (coord(mline(i,1),1)==0 | coord(mline(i,1),1)==dim_x | coord(mline(i,1),2)==0 | coord(mline(i,1),2)==dim_y)
        if (coord(mline(i,2),1)==0 | coord(mline(i,2),1)==dim_x | coord(mline(i,2),2)==0 | coord(mline(i,2),2)==dim_y)
            k=k+1;
            edge_lines(k)=i;
        else
            h=h+1;
            int_lines(h)=i;
        end    
    else
        h=h+1;
        int_lines(h)=i;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% OUTPUT 

fid=fopen('input_Gmsh.geo','w');
fprintf(fid,'%s \n','// INPUT DATA FOR GMSH');
fprintf(fid,'\n');
fprintf(fid,'lc=%i;\n',lc);
fprintf(fid,'\n');

for i=1:nb_nodes
    fprintf(fid,'Point(%i) = {%i,%i,0,lc};\n',i,coord(i,1),coord(i,2));
end
fprintf(fid,'\n');

for i=1:nline
    fprintf(fid,'Line(%i) = {%i,%i};\n',i,mline(i,1),mline(i,2));
end
fprintf(fid,'\n');

for i=1:nb_grains
    fprintf(fid,'Line Loop(%i) = {',i);
    for j=1:vertices_per_grains(i)-1
        fprintf(fid,'%i,',lineloop(i,j));
    end
    fprintf(fid,'%i};\n',lineloop(i,vertices_per_grains(i)));
end
fprintf(fid,'\n');

for i=1:nb_grains
    fprintf(fid,'Plane Surface(%i) = {%i};\n',i,i);
end
fprintf(fid,'\n');

% A tag is defined for the lines constituting the grain boundaries
fprintf(fid,'Physical Line(200) = {');
for i=1:length(int_lines)-1
    fprintf(fid,'%i,',int_lines(i));
end
fprintf(fid,'%i};\n',int_lines(length(int_lines)));
fprintf(fid,'\n');

% A tag is defined for the surfaces of the grains
fprintf(fid,'Physical Surface(300) = {');
for i=1:nb_grains-1
    fprintf(fid,'%i,',i);
end
fprintf(fid,'%i};\n',nb_grains);

if (quadrangle)
    fprintf(fid,'\n');
    fprintf(fid,'Recombine Surface {');
    for i=1:nb_grains-1
        fprintf(fid,'%i,',i);
    end
    fprintf(fid,'%i};\n',nb_grains);
end

fclose(fid);