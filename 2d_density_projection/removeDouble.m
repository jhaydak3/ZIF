% function Cord = removeDouble(Cord,plane_flag)
% % removes doubles - for example, if you are projecting onto XY plane, 
% % then this function will remove Zn atoms of the form (x1,y1,z) (x2,y2,z)
% % so that only one remains
% CordNoDoubleXY = [];
% [rows,numAtom] = size(Cord)
% for ndx = 1:numAtom
%     alreadyIn=0;
%     [rows,cols] = size(CordNoDoubleXY);
%     for ndx2 = 1:cols
%         if plane_flag == 'xy'
%             if CordNoDoubleXY(1,ndx2) == Cord(1,ndx) && CordNoDoubleXY(2,ndx2) == Cord(2,ndx)
%                 alreadyIn=1;
%                 break
%             end
%         elseif plane_flag == 'xz'
%             if CordNoDoubleXY(1,ndx2) == Cord(1,ndx) && CordNoDoubleXY(3,ndx2) == Cord(3,ndx)
%                 alreadyIn=1;
%                 break
%             end
%         elseif plane_flag == 'yz'
%             if CordNoDoubleXY(2,ndx2) == Cord(2,ndx) && CordNoDoubleXY(3,ndx2) == Cord(3,ndx)
%                 alreadyIn=1;
%                 break
%             end
%         else
%             error('Plane flag not recognized')
%         end
%     end
%     if alreadyIn == 0
%        CordNoDoubleXY = [CordNoDoubleXY Cord(:,ndx)];
%     end
% end
% 
% Cord = CordNoDoubleXY
% end
function Cord = removeDouble(Cord,plane_flag)
% removes doubles - for example, 
% then this function will remove Zn atoms of the form (x,y,z) (x,y,z)
% so that only one remains
CordNoDoubleXY = [];
[rows,numAtom] = size(Cord)
for ndx = 1:numAtom
    alreadyIn=0;
    [rows,cols] = size(CordNoDoubleXY);
    for ndx2 = 1:cols
        if CordNoDoubleXY(1,ndx2) == Cord(1,ndx) && CordNoDoubleXY(2,ndx2) == Cord(2,ndx) && CordNoDoubleXY(3,ndx2) == Cord(3,ndx)
            alreadyIn=1;
            break
        end
    end
    if alreadyIn == 0
        CordNoDoubleXY = [CordNoDoubleXY Cord(:,ndx)];
    end
end

Cord = CordNoDoubleXY
end