writerObj = VideoWriter('p5.avi');
writerObj.FrameRate=5;
open(writerObj);
for K = 1 : 98
  filename = sprintf('p5/%d.png', K);
  thisimage = imread(filename);
  writeVideo(writerObj, thisimage);
end
close(writerObj);