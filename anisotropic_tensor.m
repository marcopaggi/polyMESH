% anisotropic_tensor.m
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: MATLAB script to determine the 3D elastic constitutive tensor 
%          and the 2D plane stress counterpart  
%          for anisotropic Silicon, for any grain orientation. 
%
% Input:   (1) Miller indices of the plane [abc] with normal vector <abc>, 
%              defining the grain orientation
%  
%          (2) Elastic coefficients c1111, c1122, c2323 of the constitutive matrix 
%              C of the anisotropic grain 
%              in the standard reference system defined by axes <100>, <010>, <001> 
%          
% Output:  Stiffness matrix of the rotated polycrystal, C_3D, 
%          and its plane stress counterpart, C_2D 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                INPUT SECTION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=2;
b=1;
c=1;

c1111=165.7; % GPa
c1122=63.9;  % GPa
c2323=79.6;  % GPa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

norm=sqrt(a^2+b^2+c^2);
a=a/norm;
b=b/norm;
c=c/norm;

%  Computation of the angles thetax, thetay, thetaz
%
%  thetax : angle between [abc] and [100] planes
%  thetay : angle between [abc] and [010] planes
%  thetaz : angle between [abc] and [001] planes  

thetax=acos(a)
thetay=acos(b)
thetaz=acos(c)

% Non-zero componenents of the polycrystalline Silicon constitutive tensor 
% in the standard [100], [010], [001] frame

C(1,1,1,1)=c1111;
C(2,2,2,2)=c1111;
C(3,3,3,3)=c1111;
C(1,1,2,2)=c1122;
C(1,1,3,3)=c1122;
C(2,2,1,1)=c1122;
C(2,2,3,3)=c1122;
C(3,3,1,1)=c1122;
C(3,3,2,2)=c1122;
C(1,3,1,3)=c2323;
C(3,1,1,3)=c2323;
C(1,3,3,1)=c2323;
C(3,1,3,1)=c2323;
C(2,3,2,3)=c2323;
C(3,2,2,3)=c2323;
C(2,3,3,2)=c2323;
C(3,2,3,2)=c2323;
C(1,2,1,2)=c2323;
C(2,1,1,2)=c2323;
C(1,2,2,1)=c2323;
C(2,1,2,1)=c2323;

% Rotation matrices

R_100=[1, 0, 0;
       0, cos(thetax), -sin(thetax);
       0, sin(thetax), cos(thetax)];

R_010=[cos(thetay), 0, sin(thetay);
       0,1,0;
       -sin(thetay), 0, cos(thetay)];

R_001=[cos(thetaz), -sin(thetaz), 0;
       sin(thetaz), cos(thetaz), 0;
       0, 0, 1];
    
% Rotation matrix R
    
R=R_001*R_010*R_100;
        
% Stiffness tensor Cglobal of the rotated grain

Cglobal=zeros(3,3,3,3);
for i=1:3,
    for j=1:3,
        for k=1:3,
            for l=1:3,
                for p=1:3,
                    for q=1:3,
                        for r=1:3,
                            for s=1:3,
       Cglobal(i,j,k,l)=Cglobal(i,j,k,l)+R(p,i)*R(q,j)*R(r,k)*R(s,l)*C(p,q,r,s);     
                            end
                         end
                     end
                end
            end
        end
    end
end
    
C_3D=[
        [Cglobal(1,1,1,1) Cglobal(1,1,2,2) Cglobal(1,1,3,3) Cglobal(1,1,2,3) Cglobal(1,1,1,3) Cglobal(1,1,1,2)];
        [Cglobal(2,2,1,1) Cglobal(2,2,2,2) Cglobal(2,2,3,3) Cglobal(2,2,2,3) Cglobal(2,2,1,3) Cglobal(2,2,1,2)];
        [Cglobal(3,3,1,1) Cglobal(3,3,2,2) Cglobal(3,3,3,3) Cglobal(3,3,2,3) Cglobal(3,3,1,3) Cglobal(3,3,1,2)];
        [Cglobal(2,3,1,1) Cglobal(2,3,2,2) Cglobal(2,3,3,3) Cglobal(2,3,2,3) Cglobal(2,3,1,3) Cglobal(2,3,1,2)];
        [Cglobal(1,3,1,1) Cglobal(1,3,2,2) Cglobal(1,3,3,3) Cglobal(1,3,2,3) Cglobal(1,3,1,3) Cglobal(1,3,1,2)];
        [Cglobal(1,2,1,1) Cglobal(1,2,2,2) Cglobal(1,2,3,3) Cglobal(1,2,2,3) Cglobal(1,2,1,3) Cglobal(1,2,1,2)];
];
    
D=inv(C_3D);

D_2D=[[D(1,1) D(1,2) D(1,6)] ; [D(2,1) D(2,2) D(2,6)] ; [D(6,1) D(6,2) D(6,6)]]

C_2D=inv(D_2D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                OUTPUT SECTION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C_3D       % 3D constitutive elastic tensor for the anisotropic grain 
C_2D       % 2D plane stress constitutive elastic tensor for the anisotropic grain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
eig(C_2D)
