function [] = Aircraft_PlotterV(Aircraft, theta_z, theta_x)
figure 
hold
axis equal; 
grid
view(theta_z,theta_x)
LineWidth = 0.1;
if isfield(Aircraft,'Mass_prop') & isfield(Aircraft.Mass_prop,'CG_pos')
    plot3(-Aircraft.Mass_prop.CG_pos(1), 0, 0, 'pr','MarkerSize',8,'MarkerFaceColor','r')
    if isfield(Aircraft, 'x_mac')
        plot3(-Aircraft.x_mac, 0, 0, 'pb','MarkerSize',8,'MarkerFaceColor','none', 'MarkerEdgeColor','k', 'LineWidth', 1)
        legend('CG','AC')
    else
        legend('CG')
    end
end

xlabel('meters'); ylabel('meters');

% legend('CG')

xyzG = Aircraft.xyzG;
% text(-1.5,-0.7,0,'b_o','Fontsize',20)
% quiver3(-1.5,0,0,1,0,0,'linewidth',1.5,'Markersize',20)
% quiver3(-1.5,0,0,0,1,0,'linewidth',1.5)
% quiver3(-1.5,0,0,0,0,1,'linewidth',1.5)
% plot3(-1.5,0,0,'.','Markersize',20)

for i = 1:size(xyzG.crnr1,1)
        
    x = [xyzG.crnr1(i,1) 
        xyzG.crnr2(i,1)
        xyzG.crnr4(i,1) 
        xyzG.crnr3(i,1)
        xyzG.crnr1(i,1)];
        
        y = [xyzG.crnr1(i,2) 
        xyzG.crnr2(i,2)
        xyzG.crnr4(i,2) 
        xyzG.crnr3(i,2)
        xyzG.crnr1(i,2)];
        
        z = [xyzG.crnr1(i,3) 
        xyzG.crnr2(i,3)
        xyzG.crnr4(i,3) 
        xyzG.crnr3(i,3)
        xyzG.crnr1(i,3)];
    
    zo = -4*ones(5,1);
    yo = 4*ones(5,1);
        
    plot3(x,y,z,'LineWidth',LineWidth);
%     quiver3(xyzG.col(i,1), xyzG.col(i,2),xyzG.col(i,3), xyzG.Normal_col(i,1),xyzG.Normal_col(i,2), xyzG.Normal_col(i,3));
   
    
%     % Vortex slings figure
%     x = [10
%         xyzG.HS1(i,1) 
%         xyzG.HS2(i,1) 
%         10];
%     
%     y = [xyzG.HS1(i,2) 
%         xyzG.HS1(i,2) 
%         xyzG.HS2(i,2) 
%         xyzG.HS2(i,2) ];
%     
%     z = [xyzG.HS1(i,3) 
%         xyzG.HS1(i,3) 
%         xyzG.HS2(i,3) 
%         xyzG.HS2(i,3)];
%     plot3(x,y,z,':')
%     plot3(xyzG.col(i,1),xyzG.col(i,2),xyzG.col(i,3),'*','markersize',10)
%     if i == 1
%     text(xyzG.col(i,1)-2,xyzG.col(i,2)+5,xyzG.col(i,3),'b_k','Fontsize',25)
%     end
%     
%     xyz_bound = (xyzG.HS1(i,:)+ xyzG.HS2(i,:))*0.5;
%     plot3(xyz_bound(1),xyz_bound(2), xyz_bound(3),'x','markersize',15)
%     if i == 1
%         text(xyz_bound(1)-2,xyz_bound(2),xyz_bound(3)+1,'b_j','Fontsize',25)
%     end
        
    
    
end






        
        
        