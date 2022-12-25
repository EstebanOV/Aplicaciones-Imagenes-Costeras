


v = [0 0; 1 0; 1 2; 0 2];
f = [1 2 3 4];
figure(1)
hold on
patch('Faces',f,'Vertices',v,'FaceColor','black', 'FaceAlpha', 0.2, 'EdgeAlpha',0.2)
xlim([0 10])
ylim([0 2])
x=[0:10]
y=2*rand(11,1);
plot(x,y)