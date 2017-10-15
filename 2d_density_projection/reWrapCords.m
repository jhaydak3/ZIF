function [ nx,ny,nz, avec,bvec,cvec ] = reWrapCords(xv,yv,zv)
% Converts trajectories back to periodic boundary conditions
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

cart_cord = [xv'; yv'; zv'];
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
for ndx1 = 1:rows
   for ndx2 = 1:cols
        thisFrac= frac_cord(ndx1,ndx2);
        if thisFrac < 0
            while thisFrac < 0 
               thisFrac = thisFrac + 1 ;
            end
        elseif thisFrac >= 1
            while thisFrac >=1
               thisFrac = thisFrac - 1 ;
            end
        end
        frac_cord(ndx1,ndx2)=thisFrac;
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

