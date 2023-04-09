function [coor, conn] = TriangleMesh_mshread(mshfile)
% [coor, conn] = TriangleMesh_mshread(mshfile), read a T3 mesh from a GMSH
% generated MSH file. Note, that only MSH version 2 is supported (i.e. you
% have to configure GMSH to output type 2 files).

%ouverture du fichier text
fmid = fopen(mshfile,'r');

% version msh 2.2 0 8

tline = fgets(fmid);
tline = fgets(fmid);%format
tline = fgets(fmid);
tline = fgets(fmid);
tline = fgets(fmid);%nombre de noeud
nn=sscanf(tline,'%f');

%initialisation
coor=[];

for i=1:(nn);
    tline = fgets(fmid); %ligne des valeurs
    c=sscanf(tline,'%f %f %f'); % extrait des donn?es type (2.000000e+002  0.000000e+000  1.256287e+006) sous forme de colonne
    coor=[coor ; c(2:3)'];%prendre les coordonn?es  xy
end

conn=[];

tline = fgets(fmid);
tline = fgets(fmid);
tline = fgets(fmid);%nombre d'element
ne=sscanf(tline,'%f');

for i=1:(ne);
    tline = fgets(fmid); %ligne des valeurs
    c=sscanf(tline,'%f %f %f'); % extrait des donn?es type (2.000000e+002  0.000000e+000  1.256287e+006) sous forme de colonne
    conn=[conn ; c([6,8,7])'];%prendre les coordonn?es  xy
end


fclose(fmid);



