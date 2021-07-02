function viewFOWT(flPath)
    [fowt, ~] = readInputFile(flPath);
       
    morisonElements = fowt.floater.morisonElements;
    for ii = 1:numel(morisonElements) % loop all the cylinders specified
        %== Characteristics of the cylinder
        D  = morisonElements(ii).diam;
        n1 = morisonElements(ii).node1;      % First node
        n2 = morisonElements(ii).node2;      % Second node


        %== Parameters of the discretization (visualization purposes only)
        d_theta = pi/20; % angular discretization of the cylinders sections


        %== The circunference is first written in a coordinate system attached 
        %== to the cylinder. This coordinate system is composed of an unitary vector
        %== parallel to the cylinder axis, called t, and two vectors perpendicular
        %== to the axis, called x and y. The first one is in the xz plane and the
        %== second one is given by cross(zvec,xvec)
        x1  = (D/2) * cos(0:d_theta:2*pi); % x1 and y1 are relative to the
        y1  = (D/2) * sin(0:d_theta:2*pi); % circunference around the first node

        x2  = (D/2) * cos(0:d_theta:2*pi); % x2 and y2 are relative to the
        y2  = (D/2) * sin(0:d_theta:2*pi); % circunference around the second node


        % Definition of the local coordinate system vectors
        zvec  = (n2(:)' - n1(:)') / norm(n2 - n1); % we make sure that zvec is a row array
        if isequal(zvec, [0,0,1] ) % if Zlocal == Zglobal, REFlocal = REFglobal
            xvec = [1,0,0];
            yvec = [0,1,0];
        else % otherwise    
            yvec  = cross( [0,0,1] , zvec ) / norm( cross([0,0,1],zvec ) );
            xvec  = cross(yvec, zvec) / norm(cross(yvec,zvec));
        end


        % Matrices with the points that constitute the circunferences at the
        % extremities of the cylinder (in the LOCAL coordinate system)
        P1_local    = [ x1 ; y1 ; zeros( 1, size(x1,2) ) ];  
        P2_local    = [ x2 ; y2 ; repmat( norm(n2 - n1), 1, size(x2,2) ) ];

        % Matrix for transforming the coordinates from the local coordinate system
        % to the global coordinate system
        M_local2global = [xvec', yvec' , zvec'];

        % Matrices with the points that constitute the circunferences at the
        % extremities of the cylinder (in the GLOBAL coordinate system)
        P1_global = M_local2global*P1_local;
        P2_global = M_local2global*P2_local;


        % Points organized in a way that is easy to plot.
        % A number of "slices" of cylinder are input between the two
        % extremities in order to make a difference between what is below and
        % what is above the water surface
        npoints = morisonElements(ii).numIntPoints;
        X = zeros( npoints , size(P1_global,2) );
        Y = X;
        Z = X;
        for jj = 1 : size( P1_global,2 ) % the slices are input in this loop
            X(:,jj) = linspace( P1_global(1,jj) , P2_global(1,jj) , npoints ) + n1(1); % The global coordinate of the first node is added to include the position
            Y(:,jj) = linspace( P1_global(2,jj) , P2_global(2,jj) , npoints ) + n1(2); % position of the origin of the local coordinate system, which was not
            Z(:,jj) = linspace( P1_global(3,jj) , P2_global(3,jj) , npoints ) + n1(3); % previously taken in account 
        end

        % What is below the water line will be identified with dark blue
        % and what is above with standard blue.
        % The indeces of color_matrix are used with the color map given by 
        % the variable 'colormap' defined below the figure declaration
        color_matrix       = ones(size(Z)); 
        color_matrix(Z>0)  = 2;     

        %== Visualization
        if ii == 1         
            figure('color', 'w')
            figHandle = gca;            
        end        
        surf(figHandle, X,Y,Z, color_matrix) % lateral surface
        hold(figHandle, 'on')
        patch(figHandle, X(1,:), Y(1,:), Z(1,:), color_matrix(1,:)); % bottom surface
        hold(figHandle, 'on')
        patch(figHandle, X(end,:), Y(end,:), Z(end,:), color_matrix(end,:)); % top surface
        hold(figHandle, 'on')

        % Define a color map with (not so)dark blue in the first argument
        % and standard blue in the second
        colormap([0 0 0.6 ; 0 0 1]);    

    %     surf(X,Y,Z, color_matrix,'EdgeColor',[0 0 1],'LineStyle','-') % customize the cylinder appearance  

        % plot the nodes of the cylinder
        hold(figHandle, 'on')
        plot3(figHandle, n1(1),n1(2),n1(3), 'or', 'MarkerSize', 5, 'MarkerFace', 'r');
        plot3(figHandle, n2(1),n2(2),n2(3), 'or', 'MarkerSize', 5, 'MarkerFace', 'r');    
    end
end