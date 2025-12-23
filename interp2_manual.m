function v_interp = interp2_manual(X, Y, V, xq, yq)

    % Inputs:
    % X, Y: Grid vectors (1D arrays defining the grid)
    % V:    Value matrix (size length(Y) x length(X))
    % xq, yq: Query points (scalars or vectors)

    % Prepare output
    v_interp = zeros(size(xq));    

    % Grid spacing (Assuming uniform grid for speed)
    dx = X(2) - X(1);
    dy = Y(2) - Y(1);   

    % Loop over all query points
    for k = 1:length(xq)

        x = xq(k);
        y = yq(k);
      
        % 1. Find the index of the bottom-left corner (x1, y1)
        % This is the part YALMIP struggles with (finding the bin)        
        % Calculate raw index (continuous)
        idx_x_raw = (x - X(1)) / dx + 1;
        idx_y_raw = (y - Y(1)) / dy + 1;
       
        % Floor to get integer index
        ix = floor(idx_x_raw);
        iy = floor(idx_y_raw);        

        % Boundary Checks (Clamp to grid edges)
        ix = max(1, min(length(X)-1, ix));
        iy = max(1, min(length(Y)-1, iy));
        
        % 2. Identify the 4 Neighbors
        x1 = X(ix);    x2 = X(ix+1);
        y1 = Y(iy);    y2 = Y(iy+1);      

        Q11 = V(iy, ix);     Q21 = V(iy, ix+1);
        Q12 = V(iy+1, ix);   Q22 = V(iy+1, ix+1);

        % 3. Calculate Normalized Local Coordinates (0 to 1)
        u = (x - x1) / (x2 - x1);
        v = (y - y1) / (y2 - y1);

        % 4. Compute Weighted Sum (Bilinear Formula)
        % Interpolate in X direction first
        R1 = (1-u)*Q11 + u*Q21; % Value at Bottom edge
        R2 = (1-u)*Q12 + u*Q22; % Value at Top edge

        % Interpolate in Y direction
        val = (1-v)*R1 + v*R2;
        v_interp(k) = val;
    end
end