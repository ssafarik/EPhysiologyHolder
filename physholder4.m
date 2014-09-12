function [pt,pt3] = physholder4(thetaBd)
    % Compute the coordinates of the physholder.
    y=10;       % Length of one edge of the bowl rim.
    wall=2;     % Wall height of the bowl perimeter.
    tab = 4;    % Length of the mounting tab.
    gap = 0.10; % Space between fold-tabs and body.
    hole = 0.75;

    x = y;                    % Half the base length.
    thetaA = 15 * pi / 180;     % Shallow bowl angle.
    thetaB = thetaBd * pi / 180;     % Steep bowl angle.
    
    
    n=19;
    pt  = zeros(1,2);
    pt3 = zeros(1,3);

    % Compute stuff from the args.
    %z=2.5586063; % 2.679492 %2.5586063 gives a bowl angle of 15.0000.
    %y1 = y / (tan(pi/2-thetaB) - tan(thetaA));%z*tan(theta0);
    mA = -tan(thetaA);          % slope of line A.
    mB = tan(pi/2-thetaB);      % slope of line B.
    z = (mA*mB*y) / (mA-mB);    % Height of actual peak.
    z0 = y * tan(thetaA);       % Height of imaginary peak.
    y1 = (mA*y) / (mA-mB);
    b = sqrt(y1^2+z^2);
    a = sqrt(x^2+b^2);


    % Define points in 3D.
    pt3(1,:) = [0, -y1, z];
    pt3(2,:) = [-x, -y, 0];
    pt3(3,:) = [-x, 0, 0];
    pt3(4,:) = [0, 0, 0];
    y5 = wall*tan(thetaB);
    pt3(5,:) = [0, y5, -wall];
    pt3(6,:) = [-x, y5, -wall];
    pt3(7,:) = [0, y5+tab, -wall];
    pt3(8,:) = [-x, y5+tab, -wall];
    pt3(9,:) = [-x, y5, -wall];
    pt3(10,:) = [-x, -y, -wall];
    pt3(11,:) = [-x, -y, -wall];
    y12 = gap*tan(thetaB);
    pt3(12,:) = [-x, y12, -gap];
    pt3(13,:) = [-x, y5+tab, -wall];
    pt3(14,:) = [-x, y5+tab, -wall];
    pt3(15,:) = [-x+wall, y5+tab, -wall];
    pt3(16,:) = [-x+wall, y5, -wall];
    pt3(17,:) = [-x, -y, -gap];
    pt3(18,:) = [-x+wall, -y, -gap];
    pt3(19,:) = [-x+wall, -y, -wall];
    pt3(20,:) = [0,0,0];
    pt3(21,:) = [0,0,0];
    pt3(22,:) = [0,0,0];

    [m,n] = size(pt3);
    for i=1:m
        pt3(i+m,1) = -pt3(i,1);
        pt3(i+m,2) =  pt3(i,2);
        pt3(i+m,3) =  pt3(i,3);
    end

    fprintf ('Bowl angle 12 is %f\n', atan2(pt3(1,3)-pt3(2,3), pt3(1,2)-pt3(2,2))*180/pi)
    fprintf ('Bowl angle 13 is %f\n', atan2(z0, pt3(1,1)-pt3(3,1))*180/pi)

    
    % Define points in 2D.
    origin  = [0,0];
    h = sqrt(z^2+(y-y1)^2);      % Peak to bottom edge.
    r = sqrt(h^2+x^2);


    pt(1,:)  = [origin(1)+0, origin(2)-y1];        % Peak point.
    pt(2,:)  = [pt(1,1)-x, pt(1,2)-h];     % Base/Left point.
    theta2a  = ang2(pt, 1, 2);
    theta2b  = ang3(pt3,1,2,3);
    theta2c  = theta2a+theta2b-pi/2;

    pt(3,:)  = [pt(2,1)-y*sin(theta2c), pt(2,2)+y*cos(theta2c)];   % Top/Left point.
    pt(4,:)  = intersectingcircles(pt(3,:), pt(1,:), x, b);%[-z*sin(theta3b), z*cos(theta3b)];    % Peak/Top/Left point.
    theta314   = ang3(pt3,3,1,4);
    theta134   = ang3(pt3,1,3,4);

    pt(5,:)  = [pt(4,1)-wall*cos(theta314), pt(4,2)+wall*sin(theta314)];   % tab wall point A.
    pt(6,:)  = [pt(3,1)-wall*cos(theta314), pt(3,2)+wall*sin(theta314)];   % tab wall point B.
    pt(7,:)  = [pt(5,1)-tab*cos(theta314), pt(5,2)+tab*sin(theta314)];     % tab point A.
    pt(8,:)  = [pt(6,1)-tab*cos(theta314), pt(6,2)+tab*sin(theta314)];     % tab point B.
    pt(10,:) = [pt(2,1)-wall*cos(theta2c), pt(2,2)-wall*sin(theta2c)];  % Left wall.
    y9 = (wall*tan(thetaB));
    pt(9,:)  = [pt(10,1)-(y+y9)*sin(theta2c), pt(10,2)+(y+y9)*cos(theta2c)];
    pt(11,:) = [pt(2,1), pt(2,2)-wall];                                 % Bottom wall.
    pt(12,:) = [pt(6,1)-gap*cos(theta134), pt(6,2)-gap*sin(theta134)];    % pt3a
    pt(13,:) = [pt(12,1)-tab*cos(theta314), pt(12,2)+tab*sin(theta314)];  % Bracket.
    [x,y] = intersectionof(pt, 7, 8, 9, 10);
    pt(14,:) = [x,y];
    pt(15,:) = [pt(14,1)-(wall)*cos(theta2c), pt(14,2)-(wall)*sin(theta2c)];% Bracket tab.
    pt(16,:) = [pt(9,1)-(wall)*cos(theta2c), pt(9,2)-(wall)*sin(theta2c)];  % Bracket tab.
    pt(17,:) = [pt(2,1)-gap*cos(theta2c), pt(2,2)-gap*sin(theta2c)];    % pt2a
    pt(18,:) = [pt(17,1)+(wall)*sin(theta2c), pt(17,2)-(wall)*cos(theta2c)];    % Bracket.
    pt(19,:) = [pt(10,1)+(wall)*sin(theta2c), pt(10,2)-(wall)*cos(theta2c)];    % Bracket.

    theta13 = ang2(pt, 3, 1);
    pt(20,:) = [pt(1,1)+cos(theta13)*hole/2, pt(1,2)+sin(theta13)*hole/2];
    
    theta23 = ang2(pt, 3, 2);
    pt(21,:) = [pt(20,1)-cos(theta23), pt(20,2)-sin(theta23)];
    [x,y] = intersectionof(pt, 1, 2, 20, 21);
    pt(21,:) = [x,y];
    
    pt(22,:) = [pt(7,1)+cos(pi/4+theta314), pt(7,2)-sin(pi/4+theta314)];
    
    % Mirror the points.
    [m,n] = size(pt);
    for i=1:m
        pt(i+m,1) = -pt(i,1);
        pt(i+m,2) =  pt(i,2);
    end

    % Make asymmetric adjustments for tab overlap.
    [x,y] = intersectionof(pt, 7, 8, 7+m, 8+m);
    pt(7,:) = [x,y];
    [x,y] = intersectionof(pt, 5, 6, 5+m, 6+m);
    pt(5,:) = [x,y];
    [x,y] = intersectionof(pt, 3, 4, 3+m, 4+m);
    pt(4,:) = [x,y];

    
    figure(1);
    clf;
    hold on;
    drawit2(pt, m);
    axis off


    figure(2); 
    clf; 

    subplot(2,2,1)
    hold on; 
    drawit3(pt3, m);
    view (-180,0);
    axis off;

    subplot(2,2,2)
    hold on; 
    drawit3(pt3, m);
    view (-90,0);
    axis off;

    subplot(2,2,3)
    hold on; 
    drawit3(pt3, m);
    view (-180,90);
    axis off;

    subplot(2,2,4)
    hold on; 
    drawit3(pt3, m);
    view (-105,8);
    axis off;
    drawnow;
    
end
    
function rv=ang2(pt, i, j)
    rv = atan2(pt(i,2)-pt(j,2), pt(i,1)-pt(j,1));
end

function rv=ang3(pt3, i, ij, j)
    a = pt3(i,:)-pt3(ij,:);
    b = pt3(j,:)-pt3(ij,:);
    rv = acos(dot(a,b)/norm(a)/norm(b));
end

function rv=dist(pt3, i, j)
    rv = norm(pt3(i,:)-pt3(j,:));
end

function rv = intersectingcircles(A, B, ra, rb)
    AB = B-A;
    d = norm(AB);
    a = (rb^2 - ra^2 + d^2) / (2*d);
    rv = B - (a/d)*AB + [-sqrt(rb^2-a^2)/d * AB(2), ...
                         +sqrt(rb^2-a^2)/d * AB(1)];
end

% Intersecting point of two lines, i & j.
function [x,y] = intersectionof(pt, i1, i2, j1, j2)
    xi1 = pt(i1,1);
    yi1 = pt(i1,2);

    xi2 = pt(i2,1);
    yi2 = pt(i2,2);
    
    xj1 = pt(j1,1);
    yj1 = pt(j1,2);

    xj2 = pt(j2,1);
    yj2 = pt(j2,2);
    
    x = (xi1*xj1*yi2 - xi2*xj1*yi1 - xi1*xj2*yi2 + xi2*xj2*yi1 - xi1*xj1*yj2 + xi1*xj2*yj1 + xi2*xj1*yj2 - xi2*xj2*yj1)/(xi1*yj1 - xj1*yi1 - xi1*yj2 - xi2*yj1 + xj1*yi2 + xj2*yi1 + xi2*yj2 - xj2*yi2);
    y = (xi1*yi2*yj1 - xi2*yi1*yj1 - xi1*yi2*yj2 + xi2*yi1*yj2 - xj1*yi1*yj2 + xj2*yi1*yj1 + xj1*yi2*yj2 - xj2*yi2*yj1)/(xi1*yj1 - xj1*yi1 - xi1*yj2 - xi2*yj1 + xj1*yi2 + xj2*yi1 + xi2*yj2 - xj2*yi2);
end

function drawit3(pt3, n)
    scatter3(pt3(:,1), pt3(:,2), pt3(:,3));
    for m=[0,n]
        plot3(pt3([1+m,2+m,3+m,1+m],1), pt3([1+m,2+m,3+m,1+m],2), pt3([1+m,2+m,3+m,1+m],3)); % Triangle
        plot3(pt3([1+m,3+m,4+m,1+m],1), pt3([1+m,3+m,4+m,1+m],2), pt3([1+m,3+m,4+m,1+m],3)); % Triangle
        plot3(pt3([4+m,5+m,6+m,3+m],1), pt3([4+m,5+m,6+m,3+m],2), pt3([4+m,5+m,6+m,3+m],3)); % Top wall
        plot3(pt3([5+m,7+m,8+m,6+m],1), pt3([5+m,7+m,8+m,6+m],2), pt3([5+m,7+m,8+m,6+m],3)); % Tab
        plot3(pt3([3+m,9+m,10+m,2+m],1), pt3([3+m,9+m,10+m,2+m],2), pt3([3+m,9+m,10+m,2+m],3)); % Tab
        plot3(pt3([2+m,11+m],1), pt3([2+m,11+m],2), pt3([2+m,11+m],3)); % Bottom wall
        plot3(pt3([6+m,12+m,13+m,14+m],1), pt3([6+m,12+m,13+m,14+m],2), pt3([6+m,12+m,13+m,14+m],3)); % Bracket
        plot3(pt3([9+m,14+m,15+m,16+m,9+m],1), pt3([9+m,14+m,15+m,16+m,9+m],2), pt3([9+m,14+m,15+m,16+m,9+m],3)); % Bracket
        %plot3(pt3([19+m,20+m,21+m,10+m],1), pt3([19+m,20+m,21+m,10+m],2), pt3([19+m,20+m,21+m,10+m],3)); % Bracket
        %plot3(pt3([5+m,17+m,18+m,7+m],1), pt3([5+m,17+m,18+m,7+m],2), pt3([5+m,17+m,18+m,7+m],3)); % Bracket
    end
    plot3(pt3([2,2+m],1), pt3([2,2+m],2), pt3([2,2+m],3));
    plot3(pt3([11,11+m],1), pt3([11,11+m],2), pt3([11,11+m],3));
    for k=1:n*2
        text(pt3(k,1)+rand(1)*.5,pt3(k,2)+rand(1)*.5,pt3(k,3)+rand(1)*.5,sprintf('%d',k));
    end
    axis equal
end

function drawit2(pt, n)
    scatter(pt(:,1), pt(:,2),'r');
    for m=[0,n]
        line(pt([1+m,2+m,3+m,1+m],1), pt([1+m,2+m,3+m,1+m],2)); % Triangle
        line(pt([1+m,3+m,4+m,1+m],1), pt([1+m,3+m,4+m,1+m],2)); % Triangle
        line(pt([4+m,5+m,6+m,3+m],1), pt([4+m,5+m,6+m,3+m],2)); % Top wall
        line(pt([5+m,7+m,8+m,6+m],1), pt([5+m,7+m,8+m,6+m],2)); % Tab
        line(pt([3+m,9+m,10+m,2+m],1), pt([3+m,9+m,10+m,2+m],2)); % Tab
        line(pt([2+m,11+m],1), pt([2+m,11+m],2)); % Bottom wall
        line(pt([6+m,12+m,13+m,14+m],1), pt([6+m,12+m,13+m,14+m],2)); % Bracket
        line(pt([9+m,14+m,15+m,16+m,9+m],1), pt([9+m,14+m,15+m,16+m,9+m],2)); % Bracket
        line(pt([17+m,18+m,19+m,10+m],1), pt([17+m,18+m,19+m,10+m],2)); % Bracket
        %line(pt([5+m,17+m,18+m,7+m],1), pt([5+m,17+m,18+m,7+m],2)); % Bracket
    end
    line(pt([2,2+m],1), pt([2,2+m],2));
    line(pt([11,11+m],1), pt([11,11+m],2));
    for k=1:n*2
        text(pt(k,1)+.2,pt(k,2)+.4,sprintf('%d',k));
    end
    axis equal
end

