voids = readtable('../exported_dataFrames/voids.xlsx');
gevs = readtable('../exported_dataFrames/4lac_w_voidiness.xlsx');

% Prepare for plotting
figure;

% Define color and transparency
color = [0.5, 0.5, 0.5];  % Example color (light gray)
transparency = .1;      % Example transparency (50%)

n_faces = 10; % Number of faces for each sphere
distances = voids.cmvd_Mpc;
furthest = max(distances);
for i = 1:height(voids)
    [xmesh, ymesh, zmesh] = sphere(n_faces);

    ra = voids{i, "RAdeg"};  % Right Ascension in degrees
    de = voids{i, "DEdeg"};  % Declination in degrees
    dist = voids{i, "cmvd_Mpc"};  % Distance in Mpc
    radius = voids{i, "Reff_Mpc"};  % Effective radius in Mpc
    [x_off, y_off, z_off] = sph2cart(deg2rad(ra), deg2rad(de), dist);
    
    % Scale the sphere to the correct size
    x = xmesh * radius;
    y = ymesh * radius;
    z = zmesh * radius;
    
    % Shift the sphere to the correct location
    x = x + x;
    y = y + y;
    z = z + z;


    CO(:,:,1) = ones(n_faces).*0.5; % red
    CO(:,:,2) = ones(n_faces).*(dist/furthest); % green
    CO(:,:,3) = ones(n_faces).*(dist/furthest); % blue
    hold on;
    surf(x+x_off, y+y_off, z+z_off, CO)%, 'FaceAlpha', transparency, 'EdgeAlpha',1)
end

% draw LoS to GeV Galaxies
for  i = 2:height(gevs) % Skipping first row
    ra = gevs{i, 3};  % Right Ascension in degrees
    de = gevs{i, 4};  % Declination in degrees
    dist = gevs{i, end - 1};  % Distance in Mpc
    [x, y, z] = sph2cart(deg2rad(ra), deg2rad(de), dist);
    hold on
    plot3([0, x], [0, y], [0, z],'marker','o','color',[0.9290 0.6940 0.1250],'linewidth',3)
end


% Set axis properties
xlabel('X (Mpc)');  % Label for x-axis
ylabel('Y (Mpc)');  % Label for y-axis
zlabel('Z (Mpc)');  % Label for z-axis
title('Visualization of a Public Cosmic Void Catalog');  % Title for the plot
grid on
axis equal
daspect([1 1 1])
% colorbar(gca)