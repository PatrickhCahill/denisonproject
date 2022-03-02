vidObj = VideoWriter("test.mp4");

open(vidObj);

Z = peaks;
surf(Z)
axis tight
set(gca,'nextplot','replacechildren');

for k = 1:60
   surf(sin(2*pi*k/20)*Z,Z)
   currframe = getframe(gcf);
   writeVideo(vidObj,currframe);
end
close(vidObj)