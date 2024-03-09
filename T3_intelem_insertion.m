% T3_intelem_insertion.m
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: Routine for the insertion of interface finite elements 
%          between grain boundaries (for T3 continuum finite elements)
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

clc; 
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input data
dim_x=20;                       % dimension of the specimen in the x-direction
dim_y=20;                       % dimension of the specimen in the y-direction
f_scale=0.98;                   % scale factor to rescale the grains
isotropic_polycrystal=true;     % this parameter permits to generate meshes
                                % for anisotropic polycrystals. If it is 
                                % false, each grain will have a different
                                % material index, allowing the assignement
                                % of different constitutive laws.
                                % Otherwise, a single material index is
                                % assigned to all the grains.
                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File loading

fpu=fopen('Gmsh_elements.txt','r');     % connectivity matrix of the 2D finite elements
temp=fscanf(fpu,'%f',[5+3,inf]);
fclose(fpu);
gmsh_els=temp';
nb_elements=size(gmsh_els,1);
nb_grains=length(unique(gmsh_els(:,5)));
igr=unique(gmsh_els(:,5));

clear temp;
fpu=fopen('Gmsh_nodes.txt','r');        % nodal coordinates
temp=fscanf(fpu,'%f',[4,inf]);
fclose(fpu);
gmsh_nodes=temp';

clear temp;
fpu=fopen('Gmsh_lines.txt','r');        % connectivity matrix of the 1D finite elements (grain boundaries)
temp=fscanf(fpu,'%f',[7,inf]);
fclose(fpu);
gmsh_lines=temp';
clear temp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the barycentre of the grains
         
% Matrices Xe, Ye and Ze contain the coordinates x, y and z of the nodes
% of each element. The number of rows is equal to the number of elements.
% The number of columns is 5 (i.e. the number of nodes of a quadrangular
% element + 1, because the first node is copied in the fifth position to 
% create a closed loop useful to print the mesh)

Xe=zeros(nb_elements,4); 
Ye=zeros(nb_elements,4); 
Ze=zeros(nb_elements,4); 
edge=zeros(nb_elements,4); 
z=zeros(nb_elements,3);

for e=1:nb_elements
    for n=1:3
        Xe(e,n)=gmsh_nodes(gmsh_els(e,5+n),2);
        Ye(e,n)=gmsh_nodes(gmsh_els(e,5+n),3);
        Ze(e,n)=gmsh_nodes(gmsh_els(e,5+n),4);
        z(e,n)=0;
    
        % z is the flag for the nodes belonging to grain boundaries
        % (z=1 -> node on grain boundary)
        if(find(gmsh_lines(1:length(gmsh_lines),6)-gmsh_nodes(gmsh_els(e,n+5),1)==0))
            z(e,n)=1;
        end
        if(find(gmsh_lines(1:length(gmsh_lines),7)-gmsh_nodes(gmsh_els(e,n+5),1)==0))
            z(e,n)=1;
        end
        
        % check to identify if a node belongs to one of the edges of the
        % specimen        
        edge(e,n)=0;
        if (Xe(e,n)<=0.0001)
            edge(e,n)=1;
        elseif (Ye(e,n)<=0.0001)
            edge(e,n)=2;
        elseif (abs(Xe(e,n)-dim_x)<=0.0001)
            edge(e,n)=3;
        elseif (abs(Ye(e,n)-dim_y)<=0.0001)
            edge(e,n)=4;
        end
        
        % check to identify if a node belongs to one of the corners of the
        % specimen
        if ((Xe(e,n)<=0.0001)&&(Ye(e,n)<=0.0001))
            edge(e,n)=12;
        elseif ((Ye(e,n)<=0.0001)&&(abs(Xe(e,n)-dim_x)<=0.0001))
            edge(e,n)=23;
        elseif ((abs(Xe(e,n)-dim_x)<=0.0001)&&(abs(Ye(e,n)-dim_y)<=0.0001))
            edge(e,n)=34;
        elseif ((abs(Ye(e,n)-dim_y)<=0.0001)&&(Xe(e,n)<=0.0001))
            edge(e,n)=41;
        end
    end 
end
Xe(:,4)=Xe(:,1);
Ye(:,4)=Ye(:,1);
edge(:,4)=edge(:,1);

% Computation of the barycentre of each finite element
Be=zeros(nb_elements,3);
for e=1:nb_elements
    Be(e,1)=(Xe(e,1)+Xe(e,2)+Xe(e,3))/3;
    Be(e,2)=(Ye(e,1)+Ye(e,2)+Ye(e,3))/3;
    Be(e,3)=(Ze(e,1)+Ze(e,2)+Ze(e,3))/3;
end

% Computation of the surface of each finite element
Ae=zeros(nb_elements);
for e=1:nb_elements
    Ae(e)=polyarea(Xe(e,:),Ye(e,:));
end

% Total surface of the grains
g=1;
grain_id=gmsh_els(1,5);
Ag=zeros(1,nb_grains);

% First moment of inertia of the grains
Sum_AeXe=zeros(1,nb_grains);
Sum_AeYe=zeros(1,nb_grains);

for e=1:nb_elements
    if (gmsh_els(e,5)==grain_id)
    else
        grain_id=gmsh_els(e,5);
        g=g+1;     
    end
    
    Ag(g)=Ag(g)+Ae(e);
    Sum_AeXe(g)= Sum_AeXe(g) + Ae(e)*Be(e,1);
    Sum_AeYe(g)= Sum_AeYe(g) + Ae(e)*Be(e,2);
end

% Barycentre of the grains
Bg=zeros(nb_grains,2);
for g=1:nb_grains
    Bg(g,1)= Sum_AeXe(g)/Ag(g);
    Bg(g,2)= Sum_AeYe(g)/Ag(g);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shrinkage of the grains. This is done by shrinking each finite element
% with respect to the barycentre of the grain to which it belongs

% Matrix T contains: element_id, grain_id, old and new coordinates of the
% 3 vertices of each finite element. New coordinates are the old ones
% multiplied by the scale factor.
% T=(e g x1 y1 x2 y2 x3 y3 x1c y1c x2c y2c x3c y3c)
T=zeros(nb_elements,14);
g=1;
grain_id=gmsh_els(1,5);
for e=1:nb_elements
        if (gmsh_els(e,5)==grain_id)
        else
            grain_id=gmsh_els(e,5);
            g=g+1;
        end
        
        c1=[Xe(e,1) Ye(e,1)];
        c2=[Xe(e,2) Ye(e,2)];
        c3=[Xe(e,3) Ye(e,3)];
        
        cb=[Bg(g,1) Bg(g,2)];
        
        c1c=f_scale.*(c1-cb)+cb;
        c2c=f_scale.*(c2-cb)+cb;
        c3c=f_scale.*(c3-cb)+cb;

        T(e,:)=[e grain_id c1 c2 c3 c1c c2c c3c];
end

clear c1 c2 c3 c1c c2c c3c cb
clear Ae Ag Sum_AeXe Sum_AeYe Xe Ye Ze

% Connectivity matrix of the finite elements containing the shrunk 
% coordinates
Xsh=T(1:nb_elements,[9 11 13 9]);
Ysh=T(1:nb_elements,[10 12 14 10]);

% The nodes belonging to the edges of the specimen are moved back to their 
% original position (the one occupied before the shrinkage)
for e=1:nb_elements
    for n=1:4
        if (edge(e,n)==1)
            Xsh(e,n)=0;
        end
        if (edge(e,n)==2)
            Ysh(e,n)=0;
        end
        if (edge(e,n)==3)
            Xsh(e,n)=dim_x;
        end
        if (edge(e,n)==4)
            Ysh(e,n)=dim_y;
        end
        if (edge(e,n)==12)
            Xsh(e,n)=0;
            Ysh(e,n)=0;
        end
        if (edge(e,n)==23)
            Xsh(e,n)=dim_x;
            Ysh(e,n)=0;
        end
        if (edge(e,n)==34)
            Xsh(e,n)=dim_x;
            Ysh(e,n)=dim_y;
        end
        if (edge(e,n)==41)
            Xsh(e,n)=0;
            Ysh(e,n)=dim_y;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot

figure(1)
hold on
for e=1:nb_elements
    plot(Xsh(e,1:4),Ysh(e,1:4),'-b', 'LineWidth',1)
    hold on
end

hold on
plot(gmsh_nodes(:,2),gmsh_nodes(:,3),'xk'),hold on
plot(Be(:,1),Be(:,2),'xr'), hold on
plot(Bg(:,1),Bg(:,2),'or');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The new connectivity matrix of the grains is built

% List of nodes
nodes=0;
sum=1;
inew=0;
connect_el=zeros(nb_elements,4);
for e=1:nb_elements
    connect_el(e,1)=e;
    for n=1:3
		iold=gmsh_els(e,5+n);   % id of the original node, before shrinkage
        flag=0;
        if z(e,n)==1
            for j=1:sum-1
                if ((iold==nodes(j,6)) && (gmsh_els(e,5)==nodes(j,7)))
                    inew=nodes(j,1);
                    flag=1;
                end
            end
            if flag==0
                inew=sum;
                nodes(sum,1)=inew;          % node id
                nodes(sum,2)=Xsh(e,n);      % x-coordinate
                nodes(sum,3)=Ysh(e,n);      % y-coordinate
                nodes(sum,4)=z(e,n);        % flag grain boundary
                nodes(sum,5)=edge(e,n);     % flag node on specimen edge
                nodes(sum,6)=iold;          % id of the original node, before shrinkage
                nodes(sum,7)=gmsh_els(e,5); % id of the belonging grain
                sum=sum+1;
            end
        elseif z(e,n)==0
            for j=1:sum-1
                if (iold==nodes(j,6))
                    inew=nodes(j,1);
                    flag=1;
                end
            end
            if flag==0
            inew=sum;
            nodes(sum,1)=inew;
            nodes(sum,2)=Xsh(e,n);
            nodes(sum,3)=Ysh(e,n);
            nodes(sum,4)=z(e,n);
            nodes(sum,5)=edge(e,n);
            nodes(sum,6)=iold;
            nodes(sum,7)=gmsh_els(e,5);
            sum=sum+1;
            end
        end
        connect_el(e,n+1)=inew;
    end
end		

nb_nodes=length(nodes);

% Connectivity matrix of the grains
connect_gr=zeros(nb_grains,500);
% l12_gr contains the number of vertices of each grain
l12_gr=zeros(nb_grains,1);
centre=zeros(nb_grains,2);
for g=1:nb_grains
    teta=10*ones(nb_nodes,1);
    centre(g,1)=Bg(g,1);
    centre(g,2)=Bg(g,2);
    for j=1:length(nodes)
        if (nodes(j,7)==igr(g))
            if (nodes(j,4)==1 || nodes(j,5)~=0)
                l12_gr(g)=l12_gr(g)+1;
                deltay=nodes(j,3)-centre(g,2);
                deltax=nodes(j,2)-centre(g,1);
                tang(l12_gr(g))=deltay/deltax;
                if (deltay>0 && deltax>0)
                    teta(l12_gr(g))=atan(tang(l12_gr(g)));
                elseif (deltay>0 && deltax<0) || (deltay<0 && deltax<0)
                    teta(l12_gr(g))=atan(tang(l12_gr(g)))+pi;
                elseif (deltay<0 && deltax>0)
                    teta(l12_gr(g))=atan(tang(l12_gr(g)))+2*pi;
                end
               
                connect_gr(g,l12_gr(g,1))=nodes(j,1);
            end
        end
    end
    prov=zeros(1,500);
    for j=1:l12_gr(g)
        [C,index]=min(teta);
        prov(j)=connect_gr(g,index);
        teta(index)=10;
    end
    connect_gr(g,:)=prov(:);
    connect_gr(g,l12_gr(g,1)+1)=connect_gr(g,1);
end

clear Be gmsh_els edge teta Xsh Ysh


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Insertion of interface elements along the grain boundaries -> Construction
% of the connectivity matrix of the interface elements

% Construction of the matrix "seg" containing the interface segments 
% (loop over the grains) seg=[x1 y1 x2 y2 segment_id]
% segment_id = i*1000+j where i=grain_id and j=segment_id of the considered 
% grain
som=0;
for i=1:nb_grains
	for j=1:l12_gr(i)   % number of segments defining the contour of the grain
		seg(j+som,1:4)=[gmsh_nodes(nodes(connect_gr(i,j),6),2),gmsh_nodes(nodes(connect_gr(i,j),6),3),gmsh_nodes(nodes(connect_gr(i,j+1),6),2),gmsh_nodes(nodes(connect_gr(i,j+1),6),3)];      %coordinate originali
		segmn(j+som,1:4)=[nodes(connect_gr(i,j),2),nodes(connect_gr(i,j),3),nodes(connect_gr(i,j+1),2),nodes(connect_gr(i,j+1),3)];         %coordinate aggiornate
		seg(j+som,5)=i*500+j;
	end
	som=som+l12_gr(i);
end		
nb_seg=size(seg,1); % total number of segments

% Connectivity matrix of the segments, seg_to_nodes=[segment_id, node1_id, node2_id]
seg_to_nodes=zeros(nb_seg,3);
for i=1:nb_seg
	seg_to_nodes(i,1)=seg(i,5);
	for j=1:nb_nodes
		if(segmn(i,1)==nodes(j,2) && segmn(i,2)==nodes(j,3))
			 seg_to_nodes(i,2)=nodes(j,1);
		end
		if(segmn(i,3)==nodes(j,2) && segmn(i,4)==nodes(j,3))
			 seg_to_nodes(i,3)=nodes(j,1);
		end	
	end
end

% slave and master segments orientation check
xg1(:,1)=Bg(:,1);
yg1(:,1)=Bg(:,2);
compx=zeros(1,nb_seg);
compy=zeros(1,nb_seg);
compxg=zeros(1,nb_seg);
compyg=zeros(1,nb_seg);
scalar=zeros(1,nb_seg);

for i=1:nb_seg
    compx(i)=nodes(seg_to_nodes(i,3),2)-nodes(seg_to_nodes(i,2),2);     % deltax21=x2-x1
    compy(i)=nodes(seg_to_nodes(i,3),3)-nodes(seg_to_nodes(i,2),3);     % deltay21=y2-y1
    compxg(i)=nodes(seg_to_nodes(i,2),2)-xg1(floor(seg(i,5)/500),1);    % deltax1g=x1-xg
    compyg(i)=nodes(seg_to_nodes(i,2),3)-yg1(floor(seg(i,5)/500),1);    % deltay1g=y1-yg
	scalar(i)=-compy(i)*compxg(i)+compx(i)*compyg(i);                   % clockwise rotation of 90�
    if(scalar(i)>0)
        tempval=seg_to_nodes(i,2);
        seg_to_nodes(i,2)=seg_to_nodes(i,3);
        seg_to_nodes(i,3)=tempval;
    end
end

% Construction of the matrix "pair", where the segments are coupled  to
% constitute the interface elements
contp=1;
pair1=[];
for i=1:nb_seg
	for j=1:nb_seg
		if((seg(j,1)==seg(i,3) && seg(j,2)==seg(i,4) && seg(j,3)==seg(i,1) && seg(j,4)==seg(i,2)) || (seg(j,1)==seg(i,1) && seg(j,2)==seg(i,2) && seg(j,3)==seg(i,3) && seg(j,4)==seg(i,4)))
			pair1(contp,1:2)=[seg(i,5),seg(j,5)];
			contp=contp+1;
		end
	end
end					

% Duplicated elements are deleted from the matrix "pair"
for i=1:length(pair1)
		for j=1:length(pair1)
				if(pair1(j,1)==pair1(i,2) && pair1(j,2)==pair1(i,1))
					pair1(j,1:2)=[0, 0];					
				end
		end
end		

pair1=sortrows(pair1);

for i=2:length(pair1)
    if(pair1(i-1,1)==0 && pair1(i,1)~=0)
		 lenp=i;
    end
end		

pair=pair1(lenp:length(pair1),:);

% The connectivity matrix of the interface elements is built
% connect_int=[interface_element_id, node1, node2, node3, node4, interface_type]
contatf=0;
dx=zeros(1,length(pair));
dy=zeros(1,length(pair));
connect_int=zeros(length(pair),5);
for i=1:length(pair)	
	connect_int(i,1)=i;
	for j=1:length(seg_to_nodes)
		if(seg_to_nodes(j,1)==pair(i,1))
			 connect_int(i,2:3)=[seg_to_nodes(j,2) seg_to_nodes(j,3)];
		end
		if(seg_to_nodes(j,1)==pair(i,2))
            connect_int(i,4:5)=[seg_to_nodes(j,2) seg_to_nodes(j,3)];
		end	
	end	  	
		
	dx(i)=nodes(connect_int(i,2),2)+nodes(connect_int(i,5),2)-(nodes(connect_int(i,3),2)+nodes(connect_int(i,4),2));
	dy(i)=nodes(connect_int(i,2),3)+nodes(connect_int(i,5),3)-(nodes(connect_int(i,3),3)+nodes(connect_int(i,4),3));

    if(abs(dx(i))<1e-12 && abs(dy(i))<1e-12)
        contatf=contatf+1;
        vectf(contatf)=i+1;
        figure(1)
        hold on
        plot(nodes(connect_int(i,2:5),2),nodes(connect_int(i,2:5),3),'r')
        for k=1:4
            figure(1)
            text(nodes(connect_int(i,k+1),2),nodes(connect_int(i,k+1),3),num2str(nodes(connect_int(i,k+1),1),'%d'))
            hold on
        end    	

        old=connect_int(i,4);
        connect_int(i,4)=connect_int(i,5);
        connect_int(i,5)=old;
    end   				
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT for FEAP include files

% file NODES for FEAP
zerov(1:nb_nodes)=0;
feapnodes=[nodes(:,1),zerov',nodes(:,2),nodes(:,3),zerov'];

fid = fopen('feapnodes','w');
fprintf(fid,'coor\n');
fclose(fid);

fid = fopen('feapnodes','a');
fprintf(fid,'%g %g %12.15f %12.15f %g\n',feapnodes');
fclose(fid);

% file ELEMENTS for FEAP
zerov1(1:length(connect_el))=0;

if (isotropic_polycrystal)
    material_gr(1:length(connect_el))=1;
else
    material_gr(1:length(connect_el))=T(:,2);
end
	
feapelems=[connect_el(:,1),zerov1',material_gr',connect_el(:,2),connect_el(:,3),connect_el(:,4)];

fid = fopen('feapelems','w');
fprintf(fid,'elem\n');
fclose(fid);
fid = fopen('feapelems','a');
fprintf(fid,'%i %i %i %i %i %i\n',feapelems');
fclose(fid);

% The connectivity of the interface elements is added in the file
% "feapelems" (just after the connectivity of solid elements)

if (isotropic_polycrystal)
    material_int(1:length(connect_int))=2;
else
    material_int(1:length(connect_int))=nb_grains+1;
end

zerov2(1:length(pair))=0;
intelem=[length(connect_el)+connect_int(:,1),zerov2',material_int',connect_int(:,2),connect_int(:,3),connect_int(:,4),connect_int(:,5)];

fid = fopen('feapelems','a');
fprintf(fid,'%i %i %i %i %i %i %i\n',intelem');
fclose(fid);
