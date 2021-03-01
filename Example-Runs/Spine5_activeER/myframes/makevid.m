writerObj = VideoWriter('spine5activeER-RyR.mp4','MPEG-4');
open(writerObj);
for K = 0 : 150
  filename = sprintf('frame.%0.4d.png', K);
  thisimage = imread(filename);
  writeVideo(writerObj, thisimage);
  filename = sprintf('frame.%0.4d.png', K);
  thisimage = imread(filename);
  writeVideo(writerObj, thisimage);
end
close(writerObj);