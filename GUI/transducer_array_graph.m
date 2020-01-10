function [trans_number] = transducer_array_graph(datum, Radius)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
trans_number=length(datum);
x_coord=Radius*sin(datum/180*pi);
y_coord=Radius*cos(datum/180*pi);
% figure('color','white');
for i=1:1:(length(x_coord))
    if i==1
            plot(x_coord(i),y_coord(i),'bo','linewidth',2,'markersize',10);hold on;
    else
        
    plot(x_coord(i),y_coord(i),'ro','linewidth',2,'markersize',10);hold on;
    end
end
theta=0:0.01:2*pi;
plot(Radius*cos(theta),Radius*sin(theta),'k','linewidth',1);
for i=1:1:(length(x_coord))
    xaxis=linspace(0,x_coord(i),1000);
    yaxis=linspace(0,y_coord(i),1000);
    plot(xaxis,yaxis,':','color','k','linewidth',0.5);hold on;
end
xaxis=linspace(0,0,1000);
yaxis=linspace(-1.2*Radius,1.2*Radius,1000);
plot(xaxis,yaxis,'-.','color','b','linewidth',0.5);hold on;
xaxis=linspace(-1.2*Radius,1.2*Radius,1000);
yaxis=linspace(0,0,1000);
plot(xaxis,yaxis,'-.','color','b','linewidth',0.5);hold on;
axis off;

end

