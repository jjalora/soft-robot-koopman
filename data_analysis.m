%%% Set sampling parameters
Ts = 0.1;

[ oldfile_name , oldfile_path ] = uigetfile( 'datafiles/*.mat' ,...
    'Choose old data file...' );

[ currfile_name, currfile_path ] = uigetfile('datafiles/*.mat' ,...
    'Choose current data file...' );

old_data = load([oldfile_path, oldfile_name]);
curr_data = load([currfile_path, currfile_name]);

% Run SVD to get principle directions
y_old_centered = old_data.y - old_data.y(:, end);
y_curr_centered = curr_data.y - curr_data.y(:, end);

[~, S_old, V_old] = svds(transpose(y_old_centered), 3);
[~, S_curr, V_curr] = svds(transpose(y_curr_centered), 3);

% Visualize principle directions
customFigure; classicColors = colororder;
for ii = 1:3
plot3([0 V_old(1,ii)],[0 V_old(2,ii)],[0 V_old(3,ii)],'Linewidth',1,'Color',classicColors(ii,:))
plot3([V_old(1,ii)],[V_old(2,ii)],[V_old(3,ii)],'.','MarkerSize',12,'Color',classicColors(ii,:),'HandleVisibility','off')

plot3([0 V_curr(1,ii)],[0 V_curr(2,ii)],[0 V_curr(3,ii)],'Linewidth',1,'Color',classicColors(ii,:))
plot3([V_curr(1,ii)],[V_curr(2,ii)],[V_curr(3,ii)],'.','MarkerSize',12,'Color',classicColors(ii,:),'HandleVisibility','off')
end
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
view(3)
axis equal

% Figure out how sampling affects data characteristics
old_data.t = old_data.t';
old_data.u = old_data.u';
old_data.y = old_data.y';
old_data = data.resample(old_data, Ts);

curr_data.t = curr_data.t';
curr_data.u = curr_data.u';
curr_data.y = curr_data.y';
curr_data = data.resample(curr_data, Ts);

% Run SVD to get principle directions
old_data.y = transpose(old_data.y);
y_old_centered = old_data.y - old_data.y(:, end);

curr_data.y = transpose(curr_data.y);
y_curr_centered = curr_data.y - curr_data.y(:, end);

[~, S_old, V_old] = svds(transpose(y_old_centered), 3);
[~, S_curr, V_curr] = svds(transpose(y_curr_centered), 3);

% Visualize principle directions
customFigure; classicColors = colororder;
for ii = 1:3
plot3([0 V_old(1,ii)],[0 V_old(2,ii)],[0 V_old(3,ii)],'Linewidth',1,'Color',classicColors(ii,:))
plot3([V_old(1,ii)],[V_old(2,ii)],[V_old(3,ii)],'.','MarkerSize',12,'Color',classicColors(ii,:),'HandleVisibility','off')

plot3([0 V_curr(1,ii)],[0 V_curr(2,ii)],[0 V_curr(3,ii)],'Linewidth',1,'Color',classicColors(ii,:))
plot3([V_curr(1,ii)],[V_curr(2,ii)],[V_curr(3,ii)],'.','MarkerSize',12,'Color',classicColors(ii,:),'HandleVisibility','off')
end
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
view(3)
axis equal

% Visualize input directions
% Run SVD to get principle directions
old_data.u = transpose(old_data.u);
u_old_centered = old_data.u - old_data.u(:, end);

curr_data.u = transpose(curr_data.u);
u_curr_centered = curr_data.u - curr_data.u(:, end);

[~, ~, U_old] = svds(transpose(u_old_centered), 3);
[~, ~, U_curr] = svds(transpose(u_curr_centered), 3);

% Visualize principle directions
customFigure; classicColors = colororder;
for ii = 1:3
plot3([0 U_old(1,ii)],[0 U_old(2,ii)],[0 U_old(3,ii)],'Linewidth',1,'Color',classicColors(ii,:))
plot3([U_old(1,ii)],[U_old(2,ii)],[U_old(3,ii)],'.','MarkerSize',12,'Color',classicColors(ii,:),'HandleVisibility','off')

plot3([0 U_curr(1,ii)],[0 U_curr(2,ii)],[0 U_curr(3,ii)],'Linewidth',1,'Color',classicColors(ii,:))
plot3([U_curr(1,ii)],[U_curr(2,ii)],[U_curr(3,ii)],'.','MarkerSize',12,'Color',classicColors(ii,:),'HandleVisibility','off')
end
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
view(3)
axis equal