ZnCord = dlmread('C:\Users\jhaydak3\Google Drive\School\Sholl Group\Results\Density plots\H2\ZIF-7_DFTUNTYPED_ZnCoords.txt');
zz = dlmread('C:\Users\jhaydak3\Google Drive\School\Sholl Group\Results\Density plots\Co2\MOLECULE_msd_traj_densityxy.txt');
xx = dlmread('C:\Users\jhaydak3\Google Drive\School\Sholl Group\Results\Density plots\Co2\ZIF-7_DFTUNTYPED_xpoints.txt');
yy = dlmread('C:\Users\jhaydak3\Google Drive\School\Sholl Group\Results\Density plots\Co2\ZIF-7_DFTUNTYPED_ypoints.txt');
single_mol = dlmread('C:\Users\jhaydak3\Google Drive\School\Sholl Group\Results\Density plots\Co2\MOLECULE_msd_xyz.txt');

title_legend = 'CO_2 Trajectory XY-Plane Intensity Projection';
Mol_string = 'CO_2 Trajectory';
EXPAND_FACTOR = 1;
plane_flag = 'xy'; % xy, xz, yz
x=ZnCord(1,:);
y=ZnCord(2,:);
z=ZnCord(3,:);
hold on
if plane_flag == 'xy'
    imagesc(xx,yy,zz)
elseif plane_flag == 'xz'
    imagesc(xx,zz,yy)
elseif plane_flag == 'yz'
    imagesc(yy,zz,xx)
else
    error('Image flag not recognized')
end
colormap(flipud(hot))
sx = single_mol(:,1);
sy = single_mol(:,2);
sz = single_mol(:,3);
[rows, numZn] = size(ZnCord);





%ZnCord = removeDouble(ZnCord,plane_flag);

[znx,zny,znz] = expandCords(ZnCord(1,:), ZnCord(2,:), ZnCord(3,:), EXPAND_FACTOR,plane_flag);
if plane_flag == 'xy'
    Zn1 = plot(znx,zny,'ko')
elseif plane_flag == 'xz'
    Zn1 = plot(znx,znz,'ko')
elseif plane_flag == 'yz'
    Zn1 = plot(zny,znz,'ko')
else
    error('Image flag not recognized')
end
ZnCord = [znx; zny; znz];
ZnCord = removeDouble(ZnCord,plane_flag);

[rows,numZn] = size(ZnCord);
con_mat = zeros(numZn,numZn);

if plane_flag == 'xy'
    ndxcord1 = 1;
    ndxcord2 = 2;
    ndxcord3 = 3;
elseif plane_flag == 'xz'
    ndxcord1 = 1;
    ndxcord2 = 3;
    ndxcord3 = 2;
elseif plane_flag == 'yz'
    ndxcord1 = 2;
    ndxcord2 = 3;
    ndxcord3 = 1
else
    error('Image flag not recognized')
end
for ndx1 = 1:numZn
    big_dist = 7.1; %doesnt matter
    cutoff = .1; % for doubles with same z coordinate
    dist1 = big_dist; dist1ndx=0;
    dist2 = big_dist; dist2ndx = 0;
    dist3 = big_dist; dist3ndx = 0;
    dist4 = big_dist; dist4ndx = 0;
    for ndx2 = ndx1:numZn
        thisDist = (ZnCord(ndxcord1,ndx1) - ZnCord(ndxcord1,ndx2))^2 + (ZnCord(ndxcord2,ndx1) - ZnCord(ndxcord2,ndx2))^2 + (ZnCord(ndxcord3,ndx1) - ZnCord(ndxcord3,ndx2))^2;
        thisDist = sqrt(thisDist);
        if thisDist < dist4 && thisDist > cutoff
            dist4 = thisDist; dist4ndx = ndx2;
        end
        if dist4 < dist3
            temp = dist3;
            tempndx = dist3ndx;
            dist3 = dist4;
            dist3ndx = dist4ndx;
            dist4 = temp;
            dist4ndx = tempndx;
        end
        if dist3 < dist2
            temp = dist2;
            tempndx = dist2ndx;
            dist2 = dist3;
            dist2ndx = dist3ndx;
            dist3 = temp;
            dist3ndx = tempndx;
        end
        if dist2 < dist1
            temp = dist1;
            tempndx = dist1ndx;
            dist1 = dist2;
            dist1ndx = dist2ndx;
            dist2 = temp;
            dist2ndx = tempndx;
        end
    end
    if dist1ndx ~= 0
     con_mat(ndx1,dist1ndx)=1;
    end
    if dist2ndx ~=0
        con_mat(ndx1,dist2ndx)=1;
    end
    if dist3ndx ~= 0
        con_mat(ndx1,dist3ndx)=1;
    end
    if dist4ndx ~= 0
        con_mat(ndx1,dist4ndx)=1;
    end
    
end
colorString = 'k';
for ndx1 = 1:numZn
   for ndx2 = 1:numZn
      if con_mat(ndx1,ndx2) == 1
         x1 = ZnCord(ndxcord1,ndx1);
         x2 = ZnCord(ndxcord1,ndx2);
         y1 = ZnCord(ndxcord2,ndx1);
         y2 = ZnCord(ndxcord2,ndx2);
         L1 = plot([x1 x2],[y1 y2],colorString)
      end
   end
end

[nx,ny,nz,avec,bvec,cvec] = reWrapCords(sx,sy,sz);

start = 2500; dt = 10; stop = 4000;
sx = nx; sy = ny; sz = nz;
quiverColorString = 'm';
trajColorString = 'b.--';
if plane_flag == 'xy'
    s1 = plot(sx(start:dt:stop),sy(start:dt:stop),trajColorString)
    quiver(0,0,avec(1),avec(2),quiverColorString)
    quiver(0,0,bvec(1),bvec(2),quiverColorString)
    xlabel(['x (angstrom)'])
    ylabel(['y (angstrom)'])
elseif plane_flag == 'xz'
    s1 = plot(sx(start:dt:stop),sz(start:dt:stop),trajColorString)
    quiver(0,0,avec(1),avec(3),quiverColorString)
    quiver(0,0,cvec(1),cvec(3),quiverColorString)
    xlabel(['x (angstrom)'])
    ylabel(['z (angstrom)'])
elseif plane_flag == 'yz'
    s1 = plot(sy(start:dt:stop),sz(start:dt:stop),trajColorString)
    quiver(0,0,bvec(2),bvec(3),quiverColorString)
    quiver(0,0,cvec(2),cvec(3),quiverColorString)
    xlabel(['y (angstrom)'])
    ylabel(['z (angstrom)'])
else
    error('Image flag not recognized')
end

%legend('Zn Atom','Linker')
legend([s1, L1, Zn1] ,Mol_string,'Linker','Zinc Atom')
title(title_legend)



