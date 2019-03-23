format longEng;

global header_lines       = 6;
global star_age_col       = 3;
global log_Teff_col       = 33;
global log_L_col          = 34;
global surface_h1_col     = 55;
global surface_li_col     = 59;
global surf_avg_v_rot_col = 42;
global sz_top_radius_col  = 79;
global sz_bot_radius_col  = 80;
global sz_top_zone_col    = 83;
global sz_bot_zone_col    = 84;
global sz_top_vrot_col    = 85;
global sz_bot_vrot_col    = 86;
global mass_conv_core_col = 8;
global mb_activated_col   = 72;

global data_parent_folder = '/home/rcaballeron/MESA/workspace/sun-jupiter-system/Docs/runs/run_paper';
global filename = '1M_photosphere_history.data';
global gauss_fields = ['0g'; '3.5g'; '4g'; '4.5g'; '5g'; '5.5g'];
global rotational_vels = ['3kms'; '5kms'; '7kms'; '10kms'; '12kms'];
global colors = ['k'; 'r'; 'g'; 'b'; 'y'; 'm'; 'c'];

%Sun constants
global sun_gauss_field = 1.0; %G
global sun_age = 4.57e9; %years
global sun_A_Li7 = 1.1;
global sun_vel_rot = 2.0; %km/s



%Prepare a parsing format for reading the column positions indicated by cols
function result = get_parsing_fmt(cols)
  fmt = '';
  ini = 0;
  for i=1:length(cols)
    fmt = [fmt, repmat('%*s',1,cols(i)-ini-1),'%f'];
    ini = cols(i);
  end
  %Discard columns till the end of line
  fmt1 = '%*[^\n]';  
  result = [fmt, fmt1];
end

%Read from file filename the number of columns num_cols according to the
%line format lin_fmt and discarding the number of header lines header_lines
function result = read_matrix_from_file(filename, line_fmt, header_lines, num_cols)
  % Open file for reading
  infileID = fopen(filename, 'r');

  %Discard header lines
  i = 1;
  while i <= header_lines
    fgetl(infileID);
    i = i +1;
  end;
  
  % First read file to count number of lines with data
  no_of_lines = 0;
  while ~feof(infileID)
    no_of_lines = no_of_lines + 1;
    fgetl(infileID);
  end
  fclose(infileID);

  % Can now define arrays x and y of known length
  A = zeros(no_of_lines, num_cols);

  % Re-open the file for reading
  infileID = fopen(filename, 'r');

  %Discard header lines
  i = 1;
  while i <= header_lines
    fgetl(infileID);
    i = i +1;
  end;
  
  % Open file for reading
  % Read x and y coordinates from the file and store in arrays
  sizeA=[num_cols, no_of_lines];
  A = fscanf(infileID, line_fmt, sizeA);
  
  fclose(infileID);
  
  %Transpose the array so that A matches the orientation of the data in the file.
  result = A';

end

% Calculate the Li abundancies. The format of the matrix must be:
% A(1) = star age
% A(2) = hydrogen abundance in surface
% A(3) = lithium abundance in surface
%
% return ALi7
function result = calculate_A_Li7(A)
  ma_li7=7.016005;
  ma_h1=1.00794;

  %Function for calculating Li abundancies
  fA_Li7 = @(h1,li7) log10((li7./ma_li7)./(h1./ma_h1)) + 12;
  
  A_Li7 = zeros(length(A), 1);
  
  A_Li7(:,1) = fA_Li7(A(:,2),A(:,3));
  
  result = A_Li7;
end 


function plot_A_Li7(A, B, color, width, axis_limits)
  %Plot values
  plot(A(:,1), B(:,1), color, 'linewidth', width);

  %Axis scales
  set(gca, 'XScale', 'log');
  
  %Axis limits
  axis(axis_limits);
  
  %Axis ticks
  set(gca,'YTick',1.0:0.5:4.5);
  xticks = get (gca, "xtick"); 
  xlabels = arrayfun (@(x) sprintf ("%.2e", x), xticks, "uniformoutput", false); 
  set (gca, "xticklabel", xlabels) ;
  yticks = get (gca, "ytick"); 
  ylabels = arrayfun (@(x) sprintf ("%2.1f", x), yticks, "uniformoutput", false); 
  set (gca, "yticklabel", ylabels);

end

function plot_hr(A, color, width, ytick, axis_limits)
  %Plot values
  plot(A(:,2), A(:,3), color, 'linewidth', width);
    
  %Axis limits
  set(gca,'YTick',-1.0:ytick:3.5);
  
  %Axis ticks
  axis(axis_limits);
  xticks = get (gca, "xtick"); 
  xlabels = arrayfun (@(x) sprintf ("%1.2f", x), xticks, "uniformoutput", false); 
  set (gca, "xticklabel", xlabels) ;
  set (gca, "xdir", "reverse");
  yticks = get (gca, "ytick"); 
  ylabels = arrayfun (@(x) sprintf ("%1.2f", x), yticks, "uniformoutput", false); 
  set (gca, "yticklabel", ylabels);

end



function plot_vel_rot(A, color, width, axis_limits)
  %Plot values
  plot(A(:,1), A(:,2), color, 'linewidth', width);
  plot(A(:,1), A(:,3), color, 'linewidth', width, 'linestyle', '--');
  plot(A(:,1), A(:,4), color, 'linewidth', width, 'linestyle', ':');

  %Axis scales
  set(gca, 'XScale', 'log');  
  
  %Axis limits
  set(gca,'YTick',axis_limits(3):10:axis_limits(4));
  
  %Axis ticks
  axis(axis_limits);
  xticks = get (gca, "xtick"); 
  xlabels = arrayfun (@(x) sprintf ("%.2e", x), xticks, "uniformoutput", false); 
  set (gca, "xticklabel", xlabels) ;
  yticks = get (gca, "ytick"); 
  ylabels = arrayfun (@(x) sprintf ("%2.1f", x), yticks, "uniformoutput", false); 
  set (gca, "yticklabel", ylabels);

end

function plot_size_cz(A, color, width, ytick, axis_limits)
  %Plot values
  plot(A(:,1), A(:,2) .- A(:,3), color, 'linewidth', width);
  
  %Axis scales
  set(gca, 'XScale', 'log');  
  
  %Axis limits
  set(gca,'YTick',0:ytick:1.0);
  
  %Axis ticks
  axis(axis_limits);
  xticks = get (gca, "xtick"); 
  xlabels = arrayfun (@(x) sprintf ("%.2e", x), xticks, "uniformoutput", false); 
  set (gca, "xticklabel", xlabels) ;
  yticks = get (gca, "ytick"); 
  ylabels = arrayfun (@(x) sprintf ("%1.2f", x), yticks, "uniformoutput", false); 
  set (gca, "yticklabel", ylabels);

end


function hr_plots(gauss_fields, rotational_vels, ytick, atitle, axis_limits)
  global data_parent_folder;
  global filename;
  global star_age_col;
  global log_Teff_col;
  global log_L_col;
  global header_lines;
  global colors;
  line_width = 3;
  
  hold('on');
  labels = {};
  
  for i=1:rows(gauss_fields)
    for j=1:rows(rotational_vels)
      sub_folder = strcat(gauss_fields(i,:), '_', rotational_vels(j,:));
      full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
      
      fmt = get_parsing_fmt([star_age_col, log_Teff_col, log_L_col]);
      
      A = read_matrix_from_file(full_path, fmt, header_lines, 3);
     
      plot_hr(A, colors(mod(i*j,7),:), line_width, ytick, axis_limits);
      
      %Generate serie labels
      labels = {labels{:}, ['HR - ', gauss_fields(i,:), ',', rotational_vels(j,:)]};
    end
  end

  grid on;
  legend(labels, "location", "southeastoutside");
  legend boxoff
  xlabel('log Teff');
  ylabel('log L/log sun');
  title(atitle);

  hold('off');
  
end


  
function age_vs_li_plots(gauss_fields, rotational_vels, atitle)
  global data_parent_folder;
  global filename;
  global star_age_col;
  global surface_h1_col;
  global surface_li_col;
  global header_lines;
  global sun_age;
  global sun_A_Li7;
  global colors;
  line_width = 3;
  
    
  hold('on');
  labels = {};
  
  for i=1:rows(gauss_fields)
    for j=1:rows(rotational_vels)
      sub_folder = strcat(gauss_fields(i,:), '_', rotational_vels(j,:));
      full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
      
      fmt = get_parsing_fmt([star_age_col, surface_h1_col, surface_li_col]);
      
      A = read_matrix_from_file(full_path, fmt, header_lines, 3);
      
      B = calculate_A_Li7(A);      
      
      plot_A_Li7(A, B, colors(i*j,:), line_width, [1.0e5,1.0e10,0,4.5]);
      
      %Generate serie labels
      labels = {labels{:}, ['A(Li)=', gauss_fields(i,:), ',', rotational_vels(j,:)]};
    end
  end
  % Plot sun reference
  plot(sun_age, sun_A_Li7, '*', 'markersize', 15, 'color', [0.5,0.1,0.8]);

  grid on;
  legend(labels, "location", "southeastoutside");
  legend boxoff
  xlabel('star age');
  ylabel('A(Li7)');
  title(atitle);

  hold('off');
  
end

function age_vs_cz_size_plots(gauss_fields, rotational_vels, ytick, atitle, axis_limits)
  global data_parent_folder;
  global filename;
  global star_age_col;
  global sz_top_radius_col;
  global sz_bot_radius_col;
  global header_lines;
  global sun_age;
  global sun_vel_rot;
  global colors;
  line_width = 3;
  
  hold('on');
  labels = {};
  
  for i=1:rows(gauss_fields)
    for j=1:rows(rotational_vels)
      sub_folder = strcat(gauss_fields(i,:), '_', rotational_vels(j,:));
      full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
      
      fmt = get_parsing_fmt([star_age_col, sz_top_radius_col, sz_bot_radius_col]);
      
      A = read_matrix_from_file(full_path, fmt, header_lines, 3);
     
      plot_size_cz(A, colors(i*j,:), line_width, ytick, axis_limits);
      
      %Generate serie labels
      labels = {labels{:}, ['Conv. Zone - ', gauss_fields(i,:), ',', rotational_vels(j,:)]};
    end
  end

  grid on;
  legend(labels, "location", "southeastoutside");
  legend boxoff
  xlabel('star age (yrs)');
  ylabel('Size conv. zone (radio conv. zone/radio star)');
  title(atitle);

  hold('off');
  
end


function age_vs_vel_plots(gauss_fields, rotational_vels, atitle)
  global data_parent_folder;
  global filename;
  global star_age_col;
  global surf_avg_v_rot_col;
  global sz_top_vrot_col;
  global sz_bot_vrot_col;
  global header_lines;
  global sun_age;
  global sun_vel_rot;
  global colors;
  line_width = 3;
  
  hold('on');
  labels = {};
  
  for i=1:rows(gauss_fields)
    for j=1:rows(rotational_vels)
      sub_folder = strcat(gauss_fields(i,:), '_', rotational_vels(j,:));
      full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
      
      fmt = get_parsing_fmt([star_age_col, surf_avg_v_rot_col, sz_top_vrot_col, sz_bot_vrot_col]);
      
      A = read_matrix_from_file(full_path, fmt, header_lines, 4);
      %Calculate maximum for y axis
      %Get vel max value, divide it by 10, plus 1, multiply by 10
      ymax = (idivide(max(A(:,2)), 10, "fix") + 1) * 10;
      
      %Get vel max value, divide it by 10, minus 1, multiply by 10
      ymin = (idivide(min(A(:,2)), 10, "fix") - 1) * 10;
           
      plot_vel_rot(A, colors(i*j,:), line_width, [1.0e5, 1.0e10, ymin, ymax]);
      
      %Generate serie labels
      labels = {labels{:}, ['At surface - ', gauss_fields(i,:), ',', rotational_vels(j,:)]};
      labels = {labels{:}, ['Sup Lim CZ - ', gauss_fields(i,:), ',', rotational_vels(j,:)]};
      labels = {labels{:}, ['Inf Lim CZ - ', gauss_fields(i,:), ',', rotational_vels(j,:)]};
    end
  end
  % Plot sun reference
  plot(sun_age, sun_vel_rot, '*', 'markersize', 15, 'color', [0.5,0.1,0.8]);

  grid on;
  legend(labels, "location", "southeastoutside");
  legend boxoff
  xlabel('star age (yrs)');
  ylabel('Rotational Vel (km/s)');
  title(atitle);

  hold('off');
  
end

function plot_mb_activation(A, color, width, axis_limits)
  %Plot values
  plot(A(:,1), A(:,3), color, 'linewidth', width);
  plot(A(:,1), A(:,2), color, 'linewidth', width, 'linestyle', '--');

  %Axis scales
  set(gca, 'XScale', 'log');  
  
  %Axis limits
  set(gca,'YTick',0.0:1:axis_limits(4));
  
  %Axis ticks
  axis(axis_limits);
  xticks = get (gca, "xtick"); 
  xlabels = arrayfun (@(x) sprintf ("%.2e", x), xticks, "uniformoutput", false); 
  set (gca, "xticklabel", xlabels) ;
  yticks = get (gca, "ytick"); 
  ylabels = arrayfun (@(x) sprintf ("%2.1f", x), yticks, "uniformoutput", false); 
  set (gca, "yticklabel", ylabels);

end


function age_vs_mb_activation(gauss_fields, rotational_vels, atitle)
  global data_parent_folder;
  global filename;
  global star_age_col;
  global mb_activated_col;
  global mass_conv_core_col;
  global header_lines;
  global sun_age;
  global sun_vel_rot;
  global colors;
  line_width = 3;
  
  hold('on');
  labels = {};
  y_labels = {};
  
  %Figure size
  f = figure(1, 'position', [150, 100, 900 , 700]);
  
  for i=1:rows(gauss_fields)
    for j=1:rows(rotational_vels)
      sub_folder = strcat(gauss_fields(i,:), '_', rotational_vels(j,:));
      full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
      
      fmt = get_parsing_fmt([star_age_col, mass_conv_core_col, mb_activated_col]);
      
      A = read_matrix_from_file(full_path, fmt, header_lines, 3);
      
      %Offset the '0' value between y_offset+0 and y_offset+0.4
      y_offset = (i*j-1)*2; %An interval of 1 unit is left between each two series
      A(:,3) = (A(:,3).*0.4) .+ y_offset;
      
      %if mb_activated_col > 0 -> assign 0.6
      %if mb_activated_col > 0 -> assign 1.0
      A(find(A(:,2) == 0.0),2) = -1;
      A(find(A(:,2) > 0.0),2) = 0.6;
      A(find(A(:,2) == -1.0),2) = 1.0;
      A(:,2) = A(:,2) .+ y_offset;
      
      %Plot series
      %y_offset + 1.5 = max of coord 'y'
      %plot_mb_activation(A, colors(i*j,:), line_width, y_offset + 1.5);
      plot_mb_activation(A, colors(i*j,:), line_width, [1.0e1, 1.0e10, 0, y_offset + 1.5]);
      
      %Generate serie labels
      labels = {labels{:}, ['MB ', gauss_fields(i,:), ',', rotational_vels(j,:)]};
      labels = {labels{:}, ['CR ', gauss_fields(i,:), ',', rotational_vels(j,:)]};
      %Generate y ticks labels 0=off, 1=on
      y_labels = {y_labels{:}, 'off', 'on', 'conv.', 'rad.'};
    end
  end
  
  %Generate as many y ticks as the double of the number of series
  %y_ticks = zeros(1,2*rows(gauss_fields)*rows(rotational_vels));
  %for k=1:columns(y_ticks)
  %  y_ticks(k) = k-1;
  %end
  
  y_ticks = zeros(1,4*rows(gauss_fields)*rows(rotational_vels));
  pos = 0;
  for k=1:columns(y_ticks)
    y_ticks(k) = pos;
    
    %Calculate increments in interval 0.4, 0.2, 0.4 and 0,2
    m = mod (k,4);
    if (m == 0)
      pos = pos + 1.0; %1.0 is the empty gap between each 2 series
    elseif (m == 1)
      pos = pos + 0.4;
    elseif (m == 2)
      pos = pos + 0.2;
    elseif (m == 3)
      pos = pos + 0.4;
    end      
  end
  
  
  yticks(y_ticks);
  yticklabels(y_labels);
  
  grid on;
  legend(labels, "location", "southeastoutside");
  legend boxoff
  xlabel('star age (yrs)');
  ylabel('Magnetic braking activation');
  title(atitle);

  hold('off');
  
  H=7;W=12;
  set(f,'PaperUnits','inches');
  set(f,'PaperOrientation','portrait');
  set(f,'PaperSize',[H,W]);
  set(f,'PaperPosition',[0,0,W,H]);
  print(f,'-dpng','-color',[atitle,'.png']);
  
end



function plot_0G_var_vel()
  global gauss_fields;
  global rotational_vels; 

  age_vs_li_plots(gauss_fields(1,:), rotational_vels(2:5,:), 'A(Li7) - 0G & variable rotational velocity');
end


function plot_3_5G_var_vel()
  global gauss_fields;
  global rotational_vels; 

  age_vs_li_plots(gauss_fields(2,:), rotational_vels(3:4,:), 'A(Li7) - 3.5G & variable rotational velocity');
end

function plot_4_0G_var_vel()
  global gauss_fields;
  global rotational_vels; 

  age_vs_li_plots(gauss_fields(3,:), rotational_vels(1:5,:), 'A(Li7)- 4.0G & variable rotational velocity');
end

function plot_vel_rot_3_5G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_vel_plots(gauss_fields(2,:), rotational_vels(3:4,:), 'Rotational vel - 3.5G');
end

function plot_vel_rot_4G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_vel_plots(gauss_fields(3,:), rotational_vels(1:5,:), 'Rotational vel - 4G');
end


function plot_hr_3_5G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  hr_plots(gauss_fields(2,:), rotational_vels(1:5,:), 0.5, 'HR - 3.5G', [3.58, 3.8, -0.5, 2.2]);
end

function plot_hr_3_5G_var_vel_z_1()
  global gauss_fields;
  global rotational_vels;
  
  hr_plots(gauss_fields(2,:), rotational_vels(1:5,:), 0.1, 'HR - 3.5G', [3.75, 3.775, -0.25, 0.4]);
end


function plot_cz_size_0G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_cz_size_plots(gauss_fields(1,:), rotational_vels(1:5,:), 0.1, 'Convective zone radius - 0G', [1.0e5, 1.0e10, 0.0, 1.05]);
end

function plot_cz_size_4G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_cz_size_plots(gauss_fields(3,:), rotational_vels(1:5,:), 0.1, 'Convective zone radius - 4G', [1.0e5, 1.0e10, 0.0, 1.05]);
end


function plot_cz_size_3_5G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_cz_size_plots(gauss_fields(2,:), rotational_vels(1:5,:), 0.1, 'Convective zone radio - 3.5G', [1.0e2, 1.0e10, 0.0, 1.05]);
end

function plot_cz_size_3_5G_var_vel_z_1()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_cz_size_plots(gauss_fields(2,:), rotational_vels(1:5,:), 0.01, 'Convective zone radio - 3.5G', [1.0e7, 1.0e10, 0.25, 0.28]);
end



function plot_vel_rot_0G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_vel_plots(gauss_fields(1,:), rotational_vels(1:5,:), 'Rotational vel - 0G');
end


function plot_age_vs_mb_activation_3_5G()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_mb_activation(gauss_fields(2,:), rotational_vels(1:5,:), 'Magnetic braking activation & core radiative - 3.5G');
end


function main()
  %plot_age_vs_mb_activation_3_5G();
  %plot_0G_var_vel();  
  %plot_3_5G_var_vel();  
  %plot_4_0G_var_vel();  
  %plot_vel_rot_3_5G_var_vel();
  %plot_vel_rot_4G_var_vel();
  %plot_vel_rot_0G_var_vel();
  %plot_cz_size_3_5G_var_vel();
  %plot_cz_size_0G_var_vel();
  plot_cz_size_4G_var_vel();
  %plot_cz_size_3_5G_var_vel_z_1();
  %plot_hr_3_5G_var_vel();
  %plot_hr_3_5G_var_vel_z_1();
end


main();


 
  



% Overwrite y with new values
%filename = 'out2.data';
%outfileID = fopen(filename, 'w');
%dlmwrite('out2.data',A,'delimiter',' ','precision','%.16e');


%https://ch.mathworks.com/help/matlab/creating_plots/change-tick-marks-and-tick-labels-of-graph-1.html#d119e2453
%https://stackoverflow.com/questions/54002845/octave-add-secondary-y-axis-to-existing-plot
%https://octave.sourceforge.io/octave/function/plotyy.html
%fig = gcf();


% Open file for writing
%i = 1;
%while i <= no_of_lines
%	fprintf(outfileID, '%10.5f %10.5f', x(i), y(i));
%	i = i + 1;
%end

%fclose(outfileID);