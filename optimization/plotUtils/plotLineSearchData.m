
function [] = plotLineSearchData(xfd)

n = size(xfd,1);

figure(20)
clf;
hold on;
xmax = max(xfd(:,1))-min(xfd(:,1));
slopex = @(x) [x-xmax/10,x+xmax/10];
slopef = @(f,g) [f-g*xmax/10,f+g*xmax/10];

plot(xfd(1:n,1),xfd(1:n,2),'x')

for k = 1:n
    plot(slopex(xfd(k,1)),slopef(xfd(k,2),xfd(k,3)))
end



end

