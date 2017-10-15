function [ nx,ny,nz,avec,bvec,cvec ] = expandCords(xv,yv,zv,sf,plane_flag)
% expands coordinates
%     
a = 44.9219;
b=44.9219;
c=31.89326;
alpha = 90;
beta = 90;
gamma = 120;

avec = [a 0 0];
bvec = [b*cosd(gamma) b*sind(gamma) 0];

c3 = 1 - cosd(alpha)^2-cosd(beta)^2-cosd(gamma)^2+2*cosd(alpha)*cosd(beta)*cosd(gamma);
c3 = sqrt(c3);
c3 = c * c3 / sind(gamma);

cvec = [c*cosd(beta) c*(cosd(alpha)-cosd(gamma)*cosd(beta))/sind(gamma) c3];

cart_cord = [xv; yv; zv];
T = [ 0, 0, 0; 0, 0 ,0 ; 0, 0, 0];
T(1,1) = 1/a;
T(1,2) = -cosd(gamma) / a / sind(gamma);
omega = c3 * a * b * sind(gamma)
T(1,3) = b*c*(cosd(alpha)*cosd(gamma)-cosd(beta))/omega/sind(gamma)
T(2,2) = 1 / b/ sind(gamma);
T(2,3) = a * c * (cosd(beta)*cosd(gamma)-cosd(alpha)) / omega / sind(gamma)
T(3,3) = a * b * sind(gamma) / omega;

frac_cord = T * cart_cord;
[rows,cols] = size(frac_cord);
upper = round(sf/2);
lower = upper;
if plane_flag == 'xy'
    for ndx = 1:cols
        for ndx2 = 1
            for ndx3 = 0
                newx1 = frac_cord(1,ndx) + ndx2;
                newx2 = frac_cord(1,ndx) - ndx2;
                newy1 = frac_cord(2,ndx) + ndx3;
                newy2 = frac_cord(2,ndx) - ndx3;
                oldz = frac_cord(3,ndx);
                newvec1 = [newx1;newy1;oldz];
                newvec2 = [newx2;newy1;oldz];
                newvec3 = [newx1; newy2;oldz];
                newvec4 = [newx2; newy2;oldz];
                frac_cord = [frac_cord newvec1 newvec2 newvec3 newvec4];
            end
        end
    end
elseif plane_flag == 'xz'
    for ndx = 1:cols
        for ndx2 = 1
            for ndx3 = 0
                newx1 = frac_cord(1,ndx) + ndx2;
                newx2 = frac_cord(1,ndx) - ndx2;
                newz1 = frac_cord(3,ndx) + ndx3;
                newz2 = frac_cord(3,ndx) - ndx3;
                oldy = frac_cord(2,ndx)
                newvec1 = [newx1;oldy;newz1];
                newvec2 = [newx2;oldy;newz1];
                newvec3 = [newx1; oldy; newz2];
                newvec4 = [newx2; oldy; newz2];
                frac_cord = [frac_cord newvec1 newvec2 newvec3 newvec4];
            end
        end
    end
elseif plane_flag == 'yz'
    for ndx = 1:cols
        for ndx2 = 1
            for ndx3 = 0
                newy1 = frac_cord(2,ndx) + ndx2;
                newy2 = frac_cord(2,ndx) - ndx2;
                newz1 = frac_cord(3,ndx) + ndx3;
                newz2 = frac_cord(3,ndx) - ndx3;
                oldx = frac_cord(1,ndx) 
                newvec1 = [oldx;newy1;newz1];
                newvec2 = [oldx;newy2;newz1];
                newvec3 = [oldx; newy1; newz2];
                newvec4 = [oldx; newy2; newz2];
                frac_cord = [frac_cord newvec1 newvec2 newvec3 newvec4];
            end
        end
    end
end
%convert back
T = zeros(3,3);
T(1,1) = a;
T(1,2) = b*cosd(gamma);
T(1,3) = c*cosd(beta);
T(2,2) = b*sind(gamma);
T(2,3) = c * (cosd(alpha)-cosd(beta)*cosd(gamma))/sind(gamma);
T(3,3) = omega / a / b / sind(gamma);

ncord = T*frac_cord;
nx = ncord(1,:);
ny = ncord(2,:);
nz = ncord(3,:);
end

