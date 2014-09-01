function a = overlapDiskSegmentRectangle(diskCenter, diskRadius, diskTheta1, diskTheta2, rectangleLower, rectangleUpper)
% 
% a = overlapDiskSegmentRectangle(diskCenter, diskRadius, diskTheta1, diskTheta2, rectangleLower, rectangleUpper)
% 
% description:
%     calculates the overlap between a disk and a rectangle if center of disk lies withing the rectangle
%
% input:
%    diskCenter, diskRadius,         center and radius of disk
%    diskTheta*                      angles of segment 
%    rectangleLower, rectangleUpper  lower left and upper right corner of the rectangle


% split into two triangles 
t1 = rectangleLower; t2 = [rectangleUpper(1), rectangleLower(2)]; t3 = rectangleUpper;
a = overlapDiskSegmentTriangle(diskCenter, diskRadius, diskTheta1, diskTheta2, t1, t2, t3);

t1 = rectangleUpper; t2 = [rectangleLower(1), rectangleUpper(2)]; t3 = rectangleLower;
a = a + overlapDiskSegmentTriangle(diskCenter, diskRadius, diskTheta1, diskTheta2, t1, t2, t3);

end








