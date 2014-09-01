function a = overlapDiskRectangle(diskCenter, diskRadius, rectangleLower, rectangleUpper)
% 
% a = overlapDiskRectangle(diskCenter, diskRadius, RecatangleLower, RecatangleUpper)
% 
% description:
%     calculates the overlap between a disk and a rectangle if center of disk lies withing the rectangle
%
% input:
%    diskCenter, diskRadius,         center and radius of disk
%    rectangleLower, rectangleUpper  lower left and upper right corner of the rectangle


% split into two triangles 
t1 = rectangleLower; t2 = [rectangleUpper(1), rectangleLower(2)]; t3 = rectangleUpper;
a = overlapDiskTriangle(diskCenter, diskRadius, t1, t2, t3);

t1 = rectangleUpper; t2 = [rectangleLower(1), rectangleUpper(2)]; t3 = rectangleLower;
a = a + overlapDiskTriangle(diskCenter, diskRadius, t1, t2, t3);

end








