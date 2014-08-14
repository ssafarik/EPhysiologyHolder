function [pt,pt3] = physholder
    % Compute the coordinates of the physholder.
    z=2.5586063; % Gives a bowl angle of 15.0000.
    y=10; 
    wall=2;
    tab = 4;
    gap = 0.25;
    for theta = 10 %0:45;

        theta0 = theta * pi/180;
        n=21;
        pt  = zeros(1,2);
        pt3 = zeros(1,3);
        x = y;                    % Half the base length.

        % Compute stuff from the args.
        y1 = z*tan(theta0);
        b = sqrt(y1^2+z^2);
        a = sqrt(x^2+b^2);

        theta0y = atan(z/y);



        % Define points in 3D.
        pt3(1,:) = [0, -y1, z];
        pt3(2,:) = [-x, -y, 0];
        pt3(3,:) = [-x, 0, 0];
        pt3(4,:) = [0, 0, 0];
        y5 = wall*tan(theta0);
        pt3(5,:) = [0, y5, -wall];
        pt3(6,:) = [-x, y5, -wall];
        pt3(7,:) = [0, y5+tab, -wall];
        pt3(8,:) = [-x, y5+tab, -wall];
        pt3(9,:) = [-x, y5, -wall];
        pt3(10,:) = [-x, -y, -wall];
        pt3(11,:) = [-x, -y, -wall];
        y12 = gap*tan(theta0);
        pt3(12,:) = [-x, y12, -gap];
        pt3(13,:) = [-x, y5+tab, -wall];
        pt3(14,:) = [-x, y5+tab, -wall];
        pt3(15,:) = [-x+wall, y5+tab, -wall];
        pt3(16,:) = [-x+wall, y5, -wall];
        pt3(17,:) = [wall, y5, -wall];
        pt3(18,:) = [wall, y5+tab, -wall];
        pt3(19,:) = [-x, -y, -gap];
        pt3(20,:) = [-x+wall, -y, -gap];
        pt3(21,:) = [-x+wall, -y, -wall];

        for i=1:n
            pt3(i+n,1) = -pt3(i,1);
            pt3(i+n,2) =  pt3(i,2);
            pt3(i+n,3) =  pt3(i,3);
        end
        for i=[17,18]
            pt3(i+n,:) = [0,0,0];
        end


        % Define points in 2D.
        origin  = [0,0];
        h = sqrt(z^2+(y-y1)^2);      % Peak to bottom edge.
        r = sqrt(h^2+x^2);


        pt(1,:)  = [origin(1)+0, origin(2)-y1];        % Peak point.
        pt(2,:)  = [pt(1,1)-x, pt(1,2)-h];     % Base/Left point.
        theta2a  = ang(pt3,23,2,1);   
        theta2b  = ang(pt3,1,2,3);
        theta2c  = theta2a+theta2b-pi/2;

        pt(3,:)  = [pt(2,1)-y*sin(theta2c), pt(2,2)+y*cos(theta2c)];   % Top/Left point.
        pt(4,:)  = intersectingcircles(pt(3,:), pt(1,:), x, b);%[-z*sin(theta3b), z*cos(theta3b)];    % Peak/Top/Left point.
        theta4   = abs(atan((pt(4,1)-pt(1,1))/((pt(4,2)-pt(1,2)))));

        pt(5,:)  = [pt(4,1)-wall*sin(theta4), pt(4,2)+wall*cos(theta4)];   % tab wall point A.
        pt(6,:)  = [pt(3,1)-wall*sin(theta4), pt(3,2)+wall*cos(theta4)];   % tab wall point B.
        pt(7,:)  = [pt(5,1)-tab*sin(theta4), pt(5,2)+tab*cos(theta4)];     % tab point A.
        pt(8,:)  = [pt(6,1)-tab*sin(theta4), pt(6,2)+tab*cos(theta4)];     % tab point B.
        pt(10,:) = [pt(2,1)-wall*cos(theta2c), pt(2,2)-wall*sin(theta2c)];  % Left wall.
        y9 = (wall*tan(theta0));
        pt(9,:)  = [pt(10,1)-(y+y9)*sin(theta2c), pt(10,2)+(y+y9)*cos(theta2c)];
        pt(11,:) = [pt(2,1), pt(2,2)-wall];                                 % Bottom wall.
        pt(12,:) = [pt(3,1)-gap*cos(theta2c-theta0), pt(3,2)-gap*sin(theta2c-theta0)];    % pt3a
        pt(13,:) = [pt(12,1)-(tab+y9*cos(theta0))*sin(theta4), pt(12,2)+(tab+y9*cos(theta0))*cos(theta4)];  % Bracket.
        pt(14,:) = [pt(9,1)-(tab)*sin(theta2c), pt(9,2)+(tab)*cos(theta2c)];    % Bracket.
        pt(15,:) = [pt(14,1)-(wall)*cos(theta2c), pt(14,2)-(wall)*sin(theta2c)];% Bracket tab.
        pt(16,:) = [pt(9,1)-(wall)*cos(theta2c), pt(9,2)-(wall)*sin(theta2c)];  % Bracket tab.

        pt(17,:) = [pt(5,1)+wall*cos(theta4), pt(5,2)+wall*sin(theta4)];     % tab point A.
        pt(18,:) = [pt(7,1)+wall*cos(theta4), pt(7,2)+wall*sin(theta4)];     % tab point A.

        pt(19,:) = [pt(2,1)-gap*cos(theta2c), pt(2,2)-gap*sin(theta2c)];    % pt2a
        pt(20,:) = [pt(19,1)+(wall)*sin(theta2c), pt(19,2)-(wall)*cos(theta2c)];    % Bracket.
        pt(21,:) = [pt(10,1)+(wall)*sin(theta2c), pt(10,2)-(wall)*cos(theta2c)];    % Bracket.



        for i=1:n
            pt(i+n,1) = -pt(i,1);
            pt(i+n,2) =  pt(i,2);
        end
        for i=[17,18]
            pt(i+n,:) = [0,0];
        end


        figure(1);
        clf;
        hold on;
        scatter(pt(:,1), pt(:,2),'r');
        for m=[0,n]
            line(pt([1+m,2+m,3+m,1+m],1), pt([1+m,2+m,3+m,1+m],2)); % Triangle
            line(pt([1+m,3+m,4+m,1+m],1), pt([1+m,3+m,4+m,1+m],2)); % Triangle
            line(pt([4+m,5+m,6+m,3+m],1), pt([4+m,5+m,6+m,3+m],2)); % Top wall
            line(pt([5+m,7+m,8+m,6+m],1), pt([5+m,7+m,8+m,6+m],2)); % Tab
            line(pt([3+m,9+m,10+m,2+m],1), pt([3+m,9+m,10+m,2+m],2)); % Tab
            line(pt([2+m,11+m],1), pt([2+m,11+m],2)); % Bottom wall
            line(pt([12+m,13+m,14+m,9+m],1), pt([12+m,13+m,14+m,9+m],2)); % Bracket
            line(pt([14+m,15+m,16+m,9+m],1), pt([14+m,15+m,16+m,9+m],2)); % Bracket
            line(pt([19+m,20+m,21+m,10+m],1), pt([19+m,20+m,21+m,10+m],2)); % Bracket
            line(pt([5+m,17+m,18+m,7+m],1), pt([5+m,17+m,18+m,7+m],2)); % Bracket
        end
        line(pt([2,2+m],1), pt([2,2+m],2));
        line(pt([11,11+m],1), pt([11,11+m],2));
        for k=1:n*2
            text(pt(k,1)+.2,pt(k,2)+.4,sprintf('%d',k));
        end
        axis equal

        figure(2); clf; hold on; 
        scatter3(pt3(:,1), pt3(:,2), pt3(:,3));
        for m=[0,n]
            plot3(pt3([1+m,2+m,3+m,1+m],1), pt3([1+m,2+m,3+m,1+m],2), pt3([1+m,2+m,3+m,1+m],3)); % Triangle
            plot3(pt3([1+m,3+m,4+m,1+m],1), pt3([1+m,3+m,4+m,1+m],2), pt3([1+m,3+m,4+m,1+m],3)); % Triangle
            plot3(pt3([4+m,5+m,6+m,3+m],1), pt3([4+m,5+m,6+m,3+m],2), pt3([4+m,5+m,6+m,3+m],3)); % Top wall
            plot3(pt3([5+m,7+m,8+m,6+m],1), pt3([5+m,7+m,8+m,6+m],2), pt3([5+m,7+m,8+m,6+m],3)); % Tab
            plot3(pt3([3+m,9+m,10+m,2+m],1), pt3([3+m,9+m,10+m,2+m],2), pt3([3+m,9+m,10+m,2+m],3)); % Tab
            plot3(pt3([2+m,11+m],1), pt3([2+m,11+m],2), pt3([2+m,11+m],3)); % Bottom wall
            plot3(pt3([12+m,13+m,14+m,9+m],1), pt3([12+m,13+m,14+m,9+m],2), pt3([12+m,13+m,14+m,9+m],3)); % Bracket
            plot3(pt3([14+m,15+m,16+m,9+m],1), pt3([14+m,15+m,16+m,9+m],2), pt3([14+m,15+m,16+m,9+m],3)); % Bracket
            plot3(pt3([19+m,20+m,21+m,10+m],1), pt3([19+m,20+m,21+m,10+m],2), pt3([19+m,20+m,21+m,10+m],3)); % Bracket
            plot3(pt3([5+m,17+m,18+m,7+m],1), pt3([5+m,17+m,18+m,7+m],2), pt3([5+m,17+m,18+m,7+m],3)); % Bracket
        end
        plot3(pt3([2,2+m],1), pt3([2,2+m],2), pt3([2,2+m],3));
        plot3(pt3([11,11+m],1), pt3([11,11+m],2), pt3([11,11+m],3));
        for k=1:n*2
            %text(pt3(k,1)+rand(1)*.5,pt3(k,2)+rand(1)*.5,pt3(k,3)+rand(1)*.5,sprintf('%d',k));
        end
        axis equal
        view (-105,8);
        drawnow
    end
    fprintf ('Bowl angle is %f\n', atan2(pt3(1,3)-pt3(2,3),pt3(1,2)-pt3(2,2))*180/pi)
    %[ang(pt3,1,2,3)*180/pi, ang(pt3,2,3,1)*180/pi, dist(pt3,2,3), dist(pt3,1,3)]
    
end
    
function rv=ang(pt3, i, ij, j)
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

    