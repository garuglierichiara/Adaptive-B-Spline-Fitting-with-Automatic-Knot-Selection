function surfspline = surf_spline_load(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function surfspline = surf_spline_load(filename)
%Reads in a spline 3D surface from file
%filename --> file name with extension .db
%surfspline <-- structure of a spline surface:
%      surfspline.deguv <-- surface degree in u and in v
%      surfspline.cp <-- control point grid of size (ncpu)x(ncpv)x3
%      surfspline.ku  <-- knot vector in u
%      surfspline.kv  <-- knot vector in v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen(filename,'r');
ss=fscanf(fid,'%s',1);
ss=fscanf(fid,'%s',1);

surfspline.deguv(1)=fscanf(fid,'%5d',1);    %degree in u 
surfspline.deguv(2)=fscanf(fid,'%5d',1);    %degree in v 

ss=fscanf(fid,'%s',1);

ncpu=fscanf(fid,'%5d',1);   %num. control points in u
ncpv=fscanf(fid,'%5d',1);

ss=fscanf(fid,'%s',1);

nu=fscanf(fid,'%5d',1);    %num. knots in u
nv=fscanf(fid,'%5d',1);

surfspline.cp=zeros(ncpu,ncpv,3);
ss=fscanf(fid,'%s',1);
for i=1:ncpu
    for j=1:ncpv
        xs=fscanf(fid,'%s',1);
        ys=fscanf(fid,'%s',1);
        zs=fscanf(fid,'%s',1);
        surfspline.cp(i,j,1)=str2num(xs);
        surfspline.cp(i,j,2)=str2num(ys);
        surfspline.cp(i,j,3)=str2num(zs);
    end
end
ss=fscanf(fid,'%s',1);
for i=1:nu
    us=fscanf(fid,'%s',1);
    surfspline.ku(i)=str2num(us);
end
ss=fscanf(fid,'%s',1);
for i=1:nv
    vs=fscanf(fid,'%s',1);
    surfspline.kv(i)=str2num(vs);
end
status=fclose(fid);
