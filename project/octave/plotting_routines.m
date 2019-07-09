format longEng;
clear global;

global header_lines       = 6;

%cols position in output files
global star_age_col       = 3;
global log_abs_m_dot_col  = 7;
global mass_conv_core_col = 8;
global log_LH_col         = 17;
global log_Teff_col       = 33;
global log_L_col          = 34;
global log_R_col          = 35;
global surf_avg_v_rot_col = 42;
global center_h1_col      = 54;
global surface_h1_col     = 55;
global surface_li_col     = 59;
global mb_activated_col   = 72;
global j_dot_col          = 77;
global sz_top_radius_col  = 79;
global sz_bot_radius_col  = 80;
global sz_top_zone_col    = 83;
global sz_bot_zone_col    = 84;
global sz_top_vrot_col    = 85;
global sz_bot_vrot_col    = 86;

global data_parent_folder = '/home/rcaballeron/MESA/workspace/sun-jupiter-system/Docs/runs/run_paper';
global filename = '1M_photosphere_history.data';
global gauss_fields = ['0g'; '3.5g'; '4g'; '4.5g'; '5g'; '5.5g'];
global rotational_vels = ['0crit';'0084crit'; '014crit'; '0196crit'; '028crit'; '0336crit';];
global colors = ['k'; 'r'; 'g'; 'b'; 'y'; 'm'; 'c'];

%Sun constants
global sun_gauss_field = 1.0; %G
global sun_age = 4.57e9; %years
global sun_A_Li7 = 1.1;
global sun_vel_rot = 2.0; %km/s
global sun_T_eff = 5772; %K
global sun_L = 3.828e26; %W

%Pleiades constants
global pleiades_age = 1.0e8; %years
global pleiades_A_Li7 = 2.95;


%ZAMS
global delta_h1_center = 0.0015;
global ratio_log_LH_vs_log_L = 0.99;

%Graphic size in inches
global line_width = 3;
global W = 18;
global H = 14;
global tick_font_size = 16;
global title_font_size = 20;
global axis_font_size = 18;
global legend_font_size = 18;




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

% Calculate ZAMS age position according to the following criterium:
% 1) Stop once the central H mass fraction has been reduced by delta_h1_center 
% from its initial value.  
% 2) Move backward from this point until you find the timestep for which 
% the model first L_H/L_tot crosses ratio_log_LH_vs_log_L
%
% return age at which ZAMS is reached
function result = calculate_ZAMS(full_path)
  global star_age_col;
  global header_lines;
  global log_LH_col;
  global log_L_col;
  global center_h1_col;
  global delta_h1_center;
  global ratio_log_LH_vs_log_L;
  
  fmt = get_parsing_fmt([star_age_col, log_LH_col, log_L_col, center_h1_col]);
  A = read_matrix_from_file(full_path, fmt, header_lines, 4);
  
  %add to new columns to A initialized to 0
  %A(:,5) = diff h1_center
  %A(:,6) = log_LH div log_L
  A = [A zeros(size(A), 2)];
  
  for i=1:rows(A)
    %Calculate delta_h1
    if (i == 1)
      %For the first iteration there's no change
      A(i,5) = 0;      
    else
      %The amount of H in center decreses as time goes by
      A(i,5) = A(i-1,4) - A(i,4);
    end
    
    %Calculate ratio_log_LH_vs_log_L
    A(i,6) = power(10,A(i,2)) / power(10,A(i,3));
  end
  
  %Find those row in which diff h1 > delta_h1_center
  B = find(A(:,5) > delta_h1_center);
  %Select the first one and get rows from A till that row
  C = A(1:B(1),:);
  % Now find those row in which LH div L > ratio_log_LH_vs_log_L
  D = find(C(:,6) > ratio_log_LH_vs_log_L);
  % Select the first one
  %B(1)
  %D(1)
  %A(D(1),5)
  %A(D(1),6)
  
  %Age at which ZAMS is reached
  result = A(D(1),1);
   
end

function plot_ZAMS(A, color, width, ytick, axis_limits)
  global tick_font_size;
  
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
  set(gca, 'fontsize', tick_font_size);

end


% Calculate the Li abundancies. The format of the matrix must be:
% A(1) = star age
% A(2) = hydrogen abundance in surface
% A(3) = lithium abundance in surface
%
% return ALi7
function result = calculate_A_Li7(A)
  ma_li7 = 7.016005;
  ma_h1 = 1.00794;

  %Function for calculating Li abundancies
  fA_Li7 = @(h1,li7) log10((li7./ma_li7)./(h1./ma_h1)) + 12;
  
  A_Li7 = zeros(length(A), 1);
  
  A_Li7(:,1) = fA_Li7(A(:,2),A(:,3));
  
  result = A_Li7;
end 


function plot_A_Li7(A, B, color, width, ytick, axis_limits)
  global tick_font_size
  
  %Plot values
  plot(A(:,1), B(:,1), color, 'linewidth', width);

  %Axis scales
  set(gca, 'XScale', 'log');
  
  %Axis limits
  axis(axis_limits);
  
  %Axis ticks
  set(gca,'YTick',1.0:ytick:4.5);
  xticks = get (gca, "xtick"); 
  xlabels = arrayfun (@(x) sprintf ("%.2e", x), xticks, "uniformoutput", false); 
  set (gca, "xticklabel", xlabels) ;
  yticks = get (gca, "ytick"); 
  ylabels = arrayfun (@(x) sprintf ("%2.1f", x), yticks, "uniformoutput", false); 
  set (gca, "yticklabel", ylabels);
  set(gca, 'fontsize', tick_font_size);

end

function plot_hr(A, color, width, xy_ticks, axis_limits)
  global tick_font_size;
  %Plot values
  plot(A(:,2), A(:,3), color, 'linewidth', width);
    
  %Axis limits
  %x axis is reversed, that's why we decrement axis_limits(2)
  %in xy_ticks(1) intervals
  set(gca,'XTick',axis_limits(2):-xy_ticks(1):axis_limits(1));
  set(gca,'YTick',axis_limits(3):xy_ticks(2):axis_limits(4));  
  
  %Axis ticks
  axis(axis_limits);
  xticks = get (gca, "xtick");
  
  temp = regexp(num2str(xy_ticks(1)),'\.','split');
  precX = ['%1.' num2str(length(temp{2})) 'f'];
  
  xlabels = arrayfun (@(x) sprintf (precX, x), xticks, "uniformoutput", false); 
  set (gca, "xticklabel", xlabels) ;
  set (gca, "xdir", "reverse");
  yticks = get (gca, "ytick"); 
  ylabels = arrayfun (@(x) sprintf ("%1.2f", x), yticks, "uniformoutput", false); 
  set (gca, "yticklabel", ylabels);
  set(gca, 'fontsize', tick_font_size);
  
  %errorbar (x, y1, err, err, err, err, "~>"
  
end



function plot_vel_rot(A, color, width, ytick, axis_limits)
  global tick_font_size;
  
  %Plot values
  plot(A(:,1), A(:,2), color, 'linewidth', width);
  plot(A(:,1), A(:,3), color, 'linewidth', width, 'linestyle', '--');
  plot(A(:,1), A(:,4), color, 'linewidth', width, 'linestyle', ':');

  %Axis scales
  set(gca, 'XScale', 'log');  
  
  %Axis limits
  set(gca,'YTick',axis_limits(3):ytick:axis_limits(4));
  
  %Axis ticks
  axis(axis_limits);
  xticks = get (gca, "xtick"); 
  xlabels = arrayfun (@(x) sprintf ("%.2e", x), xticks, "uniformoutput", false); 
  set (gca, "xticklabel", xlabels) ;
  yticks = get (gca, "ytick"); 
  ylabels = arrayfun (@(x) sprintf ("%2.1f", x), yticks, "uniformoutput", false); 
  set (gca, "yticklabel", ylabels);
  set(gca, 'fontsize', tick_font_size);

end

function plot_size_cz(A, color, width, ytick, axis_limits)
  global tick_font_size;
  
  %Plot values
  plot(A(:,1), A(:,3) .- A(:,4), color, 'linewidth', width);
  %plot(A(:,1), 10 .** A(:,2), color, 'linewidth', width, 'linestyle', '--');
  
  %Axis scales
  set(gca, 'XScale', 'log');  
  
  %Axis limits
  set(gca,'YTick',0:ytick:1.0);
  
  %Axis ticks
  axis(axis_limits);
  xticks = get (gca, "xtick"); 
  xlabels = arrayfun (@(x) sprintf ("%.2e", x), xticks, "uniformoutput", false); 
  set (gca, "xticklabel", xlabels);
  yticks = get (gca, "ytick"); 
  ylabels = arrayfun (@(x) sprintf ("%1.2f", x), yticks, "uniformoutput", false); 
  set (gca, "yticklabel", ylabels);
  set(gca, 'fontsize', tick_font_size);

end

function plot_m_dot(A, color, width, ytick, axis_limits)
  global tick_font_size;
  
  %Plot values
  plot(A(:,1), A(:,2), color, 'linewidth', width);
  
  %Axis scales
  set(gca, 'XScale', 'log');  
  
  %Axis limits
  set(gca,'YTick',axis_limits(3):ytick:axis_limits(4));
  
  %Axis ticks
  axis(axis_limits);
  xticks = get (gca, "xtick"); 
  xlabels = arrayfun (@(x) sprintf ("%.2e", x), xticks, "uniformoutput", false); 
  set (gca, "xticklabel", xlabels) ;
  yticks = get (gca, "ytick"); 
  ylabels = arrayfun (@(x) sprintf ("%1.2f", x), yticks, "uniformoutput", false); 
  set (gca, "yticklabel", ylabels);
  set(gca, 'fontsize', tick_font_size);

end



function hr_plots(gauss_fields, rotational_vels, is_var_vel, xy_ticks, axis_limits, atitle, afilename)
  global data_parent_folder;
  global filename;
  global star_age_col;
  global log_Teff_col;
  global log_L_col;
  global header_lines;
  global colors;
  global sun_T_eff;
  global title_font_size;
  global axis_font_size;
  global legend_font_size;
  global line_width;
  
  hold('on');
  labels = {};
  
  f = format_figure();
  
  for i=1:rows(gauss_fields)
    for j=1:rows(rotational_vels)
      sub_folder = strcat(gauss_fields(i,:), '_', rotational_vels(j,:));
      full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
      
      fmt = get_parsing_fmt([star_age_col, log_Teff_col, log_L_col]);
      
      A = read_matrix_from_file(full_path, fmt, header_lines, 3);
     
      plot_hr(A, colors(mod(i*j,7),:), line_width, xy_ticks, axis_limits);
      
      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['HR - ', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['HR - ', strtrim(gauss_fields(i,:))]};
       endif
    end
  end
  
  % Plot sun reference
  plot(log10(sun_T_eff), 0, '*', 'markersize', 15, 'color', [0.5,0.1,0.8]);

  grid on;
  %l = legend(labels, "location", "southoutside", "orientation", "horizontal");
  l = legend(labels, "location", "southeastoutside");
  set (l, "fontsize", legend_font_size);
  legend boxoff
  xlabel('log Teff', 'fontsize', axis_font_size);
  ylabel('log (L/L_{sun})', 'fontsize', axis_font_size);
  %ylabel('{\bf \omega} = a {\bf V}');
  title(atitle, 'fontsize', title_font_size);

  hold('off');
  
  save_figure(f, afilename);
end


  
function age_vs_li_plots(gauss_fields, rotational_vels, is_var_vel, ytick, axis_limits, atitle, afilename)
  global data_parent_folder;
  global filename;
  global star_age_col;
  global surface_h1_col;
  global surface_li_col;
  global header_lines;
  global sun_age;
  global sun_A_Li7;
  global pleiades_age;
  global pleiades_A_Li7;  
  global colors;
  global title_font_size;
  global axis_font_size;
  global legend_font_size;
  global line_width;
  
    
  hold('on');
  labels = {};
  
  f = format_figure();
  
  for i=1:rows(gauss_fields)
    for j=1:rows(rotational_vels)
      sub_folder = strcat(gauss_fields(i,:), '_', rotational_vels(j,:));
      full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
       
      fmt = get_parsing_fmt([star_age_col, surface_h1_col, surface_li_col]);
      
      A = read_matrix_from_file(full_path, fmt, header_lines, 3);
      
      B = calculate_A_Li7(A);      
      
      plot_A_Li7(A, B, colors(i*j,:), line_width, ytick, axis_limits);
      
      % Plot ZAMS reference
      zams = calculate_ZAMS(full_path);
      line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 2, "linestyle", "--", "color", colors(i*j,:))
    
      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['A(Li)-', strtrim(rotational_vels(j,:))]};
        labels = {labels{:}, ['ZAMS-', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['A(Li)-', strtrim(gauss_fields(i,:))]};
        labels = {labels{:}, ['ZAMS-', strtrim(gauss_fields(i,:))]};
      endif
      
    end
  end
  % Plot sun reference
  plot(sun_age, sun_A_Li7, '*', 'markersize', 15, 'color', [0.5,0.1,0.8]);
  plot(pleiades_age, pleiades_A_Li7, 's', 'markersize', 10, 'color', [0.5,0.1,0.8], 'markerfacecolor', [0.5,0.1,0.8]);

  grid on;
  l = legend(labels, "location", "southeastoutside");
  set (l, "fontsize", legend_font_size);
  legend boxoff
  xlabel('star age', 'fontsize', axis_font_size);
  ylabel('A(Li7)', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');
  
  save_figure(f, afilename);
  
end

function age_vs_cz_size_plots(gauss_fields, rotational_vels, is_var_vel, ytick, axis_limits, atitle, afilename)
  global data_parent_folder;
  global filename;
  global star_age_col;
  global log_R_col;
  global sz_top_radius_col;
  global sz_bot_radius_col;
  global header_lines;
  global sun_age;
  global sun_vel_rot;
  global colors;
  global title_font_size;
  global axis_font_size;
  global legend_font_size;
  global line_width;
  
  hold('on');
  labels = {};
  
  f = format_figure();
  
  for i=1:rows(gauss_fields)
    for j=1:rows(rotational_vels)
      sub_folder = strcat(gauss_fields(i,:), '_', rotational_vels(j,:));
      full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
      
      fmt = get_parsing_fmt([star_age_col, log_R_col, sz_top_radius_col, sz_bot_radius_col]);
      
      A = read_matrix_from_file(full_path, fmt, header_lines, 4);
     
      plot_size_cz(A, colors(i*j,:), line_width, ytick, axis_limits);
      
      % Plot ZAMS reference
      zams = calculate_ZAMS(full_path);
      line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 2, "linestyle", "--", "color", colors(i*j,:))    
      
      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['CZ radius-', strtrim(rotational_vels(j,:))]};
        labels = {labels{:}, ['ZAMS-', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['CZ radius-', strtrim(gauss_fields(i,:))]};
        labels = {labels{:}, ['ZAMS-', strtrim(gauss_fields(i,:))]};
      endif        
    end
  end

  grid on;
  l = legend(labels, "location", "southeastoutside");
  set (l, "fontsize", legend_font_size);
  legend boxoff
  xlabel('star age (yrs)', 'fontsize', axis_font_size);
  ylabel('Size conv. zone (radius conv. zone/radius star)', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');
  
  save_figure(f, afilename);
end

function age_vs_m_dot_plots(gauss_fields, rotational_vels, is_var_vel, ytick, axis_limits, atitle, afilename)
  global data_parent_folder;
  global filename;
  global star_age_col;
  global log_abs_m_dot_col;
  global header_lines;
  global sun_age;
  global colors;
  global title_font_size;
  global axis_font_size;
  global legend_font_size;
  global line_width;
  
  hold('on');
  labels = {};
  
  f = format_figure();
  
  for i=1:rows(gauss_fields)
    for j=1:rows(rotational_vels)
      sub_folder = strcat(gauss_fields(i,:), '_', rotational_vels(j,:));
      full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
      
      fmt = get_parsing_fmt([star_age_col, log_abs_m_dot_col]);
      
      A = read_matrix_from_file(full_path, fmt, header_lines, 2);
     
      plot_m_dot(A, colors(i*j,:), line_width, ytick, axis_limits);
      
      % Plot ZAMS reference
      zams = calculate_ZAMS(full_path);
      line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 2, "linestyle", "--", "color", colors(i*j,:))    
      
      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['Mass loss-', strtrim(rotational_vels(j,:))]};
        labels = {labels{:}, ['ZAMS-', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['Mass loss-', strtrim(gauss_fields(i,:))]};
        labels = {labels{:}, ['ZAMS-', strtrim(gauss_fields(i,:))]};
      endif
      
    end
  end

  grid on;
  l = legend(labels, "location", "southeastoutside");
  set (l, "fontsize", legend_font_size);
  legend boxoff
  xlabel('star age (yrs)', 'fontsize', axis_font_size);
  ylabel('Mass loss (log(solar mass/yr))', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');
  
  save_figure(f, afilename);  
end



function age_vs_vel_plots(gauss_fields, rotational_vels, is_var_vel, ytick, x_limits, atitle, afilename)
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
  global title_font_size;
  global axis_font_size;
  global legend_font_size;
  global line_width;
  
  hold('on');
  labels = {};
  
  f = format_figure();
 
  for i=1:rows(gauss_fields)
    for j=1:rows(rotational_vels)
      sub_folder = strcat(gauss_fields(i,:), '_', rotational_vels(j,:));
      full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
      
      fmt = get_parsing_fmt([star_age_col, surf_avg_v_rot_col, sz_top_vrot_col, sz_bot_vrot_col]);
      
      A = read_matrix_from_file(full_path, fmt, header_lines, 4);
      %Get the index of the last records lower than or equal to the temporal limits
      ix_ini = find(A(:,1)<=x_limits(1), 1, 'last');
      ix_end = find(A(:,1)<=x_limits(2), 1, 'last');
      
      %Calculate maximum for y axis
      %Get vel max value, divide it by 10, plus 1, multiply by 10
      ymax = (idivide(max(A(ix_ini:ix_end,2)), ytick, "fix") + 1) * ytick;
      
      %Get vel max value, divide it by 10, minus 1, multiply by 10
      ymin = (idivide(min(A(ix_ini:ix_end,2)), ytick, "fix") - 1) * ytick;
           
      plot_vel_rot(A, colors(i*j,:), line_width, ytick, [x_limits(1), x_limits(2), ymin, ymax]);
      
      % Plot ZAMS reference
      zams = calculate_ZAMS(full_path);
      % Here we can properly assing the ymin and ymax values for all the plots
      % Each of them will potentially have a different one. We opt for fixing
      % the limits by hand.
      line("xdata",[zams,zams], "ydata",[-10,200], "linewidth", 2, "linestyle", "--", "color", colors(i*j,:))    
      
      
      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['Surface-', strtrim(rotational_vels(j,:))]};
        labels = {labels{:}, ['Top CZ-', strtrim(rotational_vels(j,:))]};
        labels = {labels{:}, ['Bottom CZ-', strtrim(rotational_vels(j,:))]};
        labels = {labels{:}, ['ZAMS-', strtrim(rotational_vels(j,:))]};      
      else
        labels = {labels{:}, ['Surface-', strtrim(gauss_fields(i,:))]};
        labels = {labels{:}, ['Top CZ-', strtrim(gauss_fields(i,:))]};
        labels = {labels{:}, ['Bottom CZ-', strtrim(gauss_fields(i,:))]};
        labels = {labels{:}, ['ZAMS-', strtrim(gauss_fields(i,:))]};            
      endif
      
    end
  end
  % Plot sun reference
  plot(sun_age, sun_vel_rot, '*', 'markersize', 15, 'color', [0.5,0.1,0.8]);
  
  grid on;
  l = legend(labels, "location", "southeastoutside");
  %l = legend(labels, "location", "southoutside", "orientation", "horizontal");
  set (l, "fontsize", legend_font_size);
  legend boxoff
  xlabel('star age (yrs)', 'fontsize', axis_font_size);
  ylabel('Rotational Vel (km/s)', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);
  
  hold('off');  
  save_figure(f, afilename);
  
end

function result = format_figure()
  global H;
  global W;

  f = figure(1, 'position', [150, 100, 900, 700]);
  set(f,'PaperUnits','inches');
  set(f,'PaperOrientation','portrait');
  set(f,'PaperSize', [H W]);
  set(f,'PaperPosition',[0 0 W H]);
  
  result = f;
  
end

function save_figure(f, title)  
  print(f,'-deps','-color',[title,'.eps']);
  close;
  %print(f,'-dpng','-color',[title,'.png']);
end

function plot_mb_activation(A, color, width, axis_limits)
  global tick_font_size;
  
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
  set(gca, 'fontsize', tick_font_size);

end


function age_vs_mb_activation(gauss_fields, rotational_vels, is_var_vel, atitle, afilename)
  global data_parent_folder;
  global filename;
  global star_age_col;
  global mb_activated_col;
  global mass_conv_core_col;
  global header_lines;
  global sun_age;
  global sun_vel_rot;
  global colors;
  global title_font_size;
  global axis_font_size;
  global legend_font_size;
  global line_width;
  
  hold('on');
  labels = {};
  y_labels = {};
  
  %Figure size
  f = format_figure();
  
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
      
      % Plot ZAMS reference
      zams = calculate_ZAMS(full_path);
      % Here we can properly assing the ymin and ymax values for all the plots
      % Each of them will potentially have a different one. We opt for fixing
      % the limits by hand.
      line("xdata",[zams,zams], "ydata",[-10,200], "linewidth", 2, "linestyle", "--", "color", colors(i*j,:))    
      
      
      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['MB-', strtrim(rotational_vels(j,:))]};
        labels = {labels{:}, ['CR-', strtrim(rotational_vels(j,:))]};
        labels = {labels{:}, ['ZAMS-', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['MB-', strtrim(gauss_fields(i,:))]};
        labels = {labels{:}, ['CR-', strtrim(gauss_fields(i,:))]};
        labels = {labels{:}, ['ZAMS-', strtrim(gauss_fields(i,:))]};
      endif
      
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
  l = legend(labels, "location", "southeastoutside");
  set (l, "fontsize", legend_font_size);
  legend boxoff
  xlabel('star age (yrs)', 'fontsize', axis_font_size);
  ylabel('Magnetic braking activation & Radiative vs. Convective core', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');

  save_figure(f, afilename);  
end



function plot_0G_var_vel()
  global gauss_fields;
  global rotational_vels; 

  age_vs_li_plots(gauss_fields(1,:), rotational_vels(1:6,:), true, 0.5, [1.0e5,1.0e10,0,4.5], 'A(Li7) - 0G & var. rotational velocity', 'li_var_vel_0_0g');
end


function plot_0G_var_vel_z1()
  global gauss_fields;
  global rotational_vels; 

  age_vs_li_plots(gauss_fields(1,:), rotational_vels(1:6,:), true, 0.1, [1.0e7, 1.0e8, 2.0, 2.5], 'A(Li7) - 0G & var. rotational velocity', 'li_var_vel_0_0g_z1');
end


function plot_3_5G_var_vel()
  global gauss_fields;
  global rotational_vels; 

  age_vs_li_plots(gauss_fields(2,:), rotational_vels(2:6,:), true, 0.5, [1.0e5,1.0e10,0,4.5], 'A(Li7) - 3.5G & var. rotational velocity', 'li_var_vel_3_5g');
end

function plot_4_0G_var_vel()
  global gauss_fields;
  global rotational_vels; 

  age_vs_li_plots(gauss_fields(3,:), rotational_vels(2:6,:), true, 0.5, [1.0e5,1.0e10,0,4.5], 'A(Li7) - 4.0G & var. rotational velocity', 'li_var_vel_4_0g');
end

function plot_4_5G_var_vel()
  global gauss_fields;
  global rotational_vels; 

  age_vs_li_plots(gauss_fields(4,:), rotational_vels(2:6,:), true, 0.5, [1.0e5,1.0e10,0,4.5], 'A(Li7) - 4.5G & var. rotational velocity', 'li_var_vel_4_5g');
end

function plot_5_0G_var_vel()
  global gauss_fields;
  global rotational_vels; 

  age_vs_li_plots(gauss_fields(5,:), rotational_vels(2:6,:), true, 0.5, [1.0e5,1.0e10,0,4.5], 'A(Li7) - 5.0G & var. rotational velocity', 'li_var_vel_5_0g');
end

function plot_5_5G_var_vel()
  global gauss_fields;
  global rotational_vels; 

  age_vs_li_plots(gauss_fields(6,:), rotational_vels(2:6,:), true, 0.5, [1.0e5,1.0e10,0,4.5], 'A(Li7) - 5.5G & var. rotational velocity', 'li_var_vel_5_5g');
end


function plot_0084vc_var_g()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_li_plots(gauss_fields(1:5,:), rotational_vels(2,:), false, 0.5, [1.0e5,1.0e10,0,4.5], 'A(Li7) - vcrit=0.0084 & var. magnetic field', 'li_vc_0084_var_g');
end

function plot_014vc_var_g()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_li_plots(gauss_fields(1:5,:), rotational_vels(3,:), false, 0.5, [1.0e5,1.0e10,0,4.5], 'A(Li7) - vcrit=0.014 & var. magnetic field', 'li_vc_014_var_g');
end
function plot_0196vc_var_g()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_li_plots(gauss_fields(1:5,:), rotational_vels(4,:), false, 0.5, [1.0e5,1.0e10,0,4.5], 'A(Li7) - vcrit=0.0196 & var. magnetic field', 'li_vc_0196_var_g');
end
function plot_028vc_var_g()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_li_plots(gauss_fields(1:5,:), rotational_vels(5,:), false, 0.5, [1.0e5,1.0e10,0,4.5], 'A(Li7) - vcrit=0.028 & var. magnetic field', 'li_vc_028_var_g');
end

function plot_0336vc_var_g()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_li_plots(gauss_fields(1:5,:), rotational_vels(6,:), false, 0.5, [1.0e5,1.0e10,0,4.5], 'A(Li7) - vcrit=0.0336 & var. magnetic field', 'li_vc_0336_var_g');
end



function plot_vel_rot_0G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_vel_plots(gauss_fields(1,:), rotational_vels(2:6,:), true, 10, [1.0e5,1.0e10], 'Rotational vel - 0G & var. rotational velocity', 'rot_vel_var_vel_0_0g');
end



function plot_vel_rot_3_5G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_vel_plots(gauss_fields(2,:), rotational_vels(2:6,:), true, 10, [1.0e5,1.0e10], 'Rotational velocity - 3.5G & var. rotational velocity', 'rot_vel_var_vel_3_5g');
end

function plot_vel_rot_4G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_vel_plots(gauss_fields(3,:), rotational_vels(2:6,:), true, 10, [1.0e5,1.0e10], 'Rotational velocity - 4.0G & var. rotational velocity', 'rot_vel_var_vel_4_0g');
end

function plot_vel_rot_4G_var_vel_z1()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_vel_plots(gauss_fields(3,:), rotational_vels(2:6,:), true, 2, [3.0e9, 1.0e10], 'Rotational velocity - 4.0G & var. rotational velocity', 'rot_vel_var_vel_4_0g_z1');
end


function plot_vel_rot_4_5G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_vel_plots(gauss_fields(4,:), rotational_vels(2:6,:), true, 10, [1.0e5,1.0e10], 'Rotational velocity - 4.5G & var. rotational velocity', 'rot_vel_var_vel_4_5g');
end

function plot_vel_rot_5G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_vel_plots(gauss_fields(5,:), rotational_vels(2:6,:), true, 10, [1.0e5,1.0e10], 'Rotational velocity - 5.0G & var. rotational velocity', 'rot_vel_var_vel_5_0g');
end

function plot_vel_rot_5_5G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_vel_plots(gauss_fields(6,:), rotational_vels(2:6,:), true, 10, [1.0e5,1.0e10], 'Rotational velocity - 5.5G & var. rotational velocity', 'rot_vel_var_vel_5_5g');
end


function plot_vel_rot_0084vc_var_g()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_vel_plots(gauss_fields(1:5,:), rotational_vels(2,:), false, 10, [1.0e5,1.0e10], 'Rotational velocity - vcrit=0.0084 & var. magnetic field', 'rot_vel_vc_0084_var_g');
end

function plot_hr_0G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  hr_plots(gauss_fields(1,:), rotational_vels(1:6,:), true, [0.05,0.5], [3.58, 3.8, -0.5, 2.2], 'HR - 0G & var. rotational velocity', 'hr_var_vel_0_0g');
end

function plot_hr_0G_var_vel_z1()
  global gauss_fields;
  global rotational_vels;
  
  hr_plots(gauss_fields(1,:), rotational_vels(1:6,:), true, [0.005,0.1], [3.75, 3.775, -0.3, 0.45], 'HR - 0G & var. rotational velocity', 'hr_var_vel_0_0g_z1');
end

function plot_hr_3_5G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  hr_plots(gauss_fields(2,:), rotational_vels(2:6,:), true, [0.05,0.5], [3.58, 3.8, -0.5, 2.2], 'HR - 3.5G & var. rotational velocity', 'hr_var_vel_3_5g');
end

function plot_hr_3_5G_var_vel_z_1()
  global gauss_fields;
  global rotational_vels;
  
  hr_plots(gauss_fields(2,:), rotational_vels(2:6,:), true, [0.005,0.1], [3.75, 3.775, -0.25, 0.4]), 'HR - 3.5G & var. rotational velocity', 'hr_var_vel_3_5g_z1';
end

function plot_hr_5_0G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  hr_plots(gauss_fields(5,:), rotational_vels(2:6,:), true, [0.05,0.5], [3.58, 3.8, -0.5, 2.2], 'HR - 5.0G & var. rotational velocity', 'hr_var_vel_5_0g');
end

function plot_hr_0336vc_var_g()
  global gauss_fields;
  global rotational_vels;
  
  hr_plots(gauss_fields(1:6,:), rotational_vels(6,:), false, [0.05,0.5], [3.58, 3.8, -0.5, 2.2], 'HR - vcrit=0.0336 & var. magnetic field', 'hr_vc_0336_var_g');
end

function plot_hr_0336vc_var_g_z1()
  global gauss_fields;
  global rotational_vels;
  
  hr_plots(gauss_fields(1:6,:), rotational_vels(6,:), false, [0.005,0.1], [3.745, 3.775, -0.3, 0.45], 'HR - vcrit=0.0336 & var. magnetic field', 'hr_vc_0336_var_g_z1');
end



function plot_cz_size_0G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_cz_size_plots(gauss_fields(1,:), rotational_vels(2:6,:), true, 0.1, [1.0e5, 1.0e10, 0.0, 1.05], 'Convective zone radius - 0G & var. rotational velocity', 'cz_var_vel_0g');
end

function plot_cz_size_3_5G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_cz_size_plots(gauss_fields(2,:), rotational_vels(2:6,:), true, 0.1, [1.0e2, 1.0e10, 0.0, 1.05], 'Convective zone size - 3.5G & var. rotational velocity', 'cz_var_vel_3_5g');
end
function plot_cz_size_4_0G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_cz_size_plots(gauss_fields(3,:), rotational_vels(2:6,:), true, 0.1, [1.0e2, 1.0e10, 0.0, 1.05], 'Convective zone size - 4.0G & var. rotational velocity', 'cz_var_vel_4_0g');
end
function plot_cz_size_4_5G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_cz_size_plots(gauss_fields(4,:), rotational_vels(2:6,:), true, 0.1, [1.0e2, 1.0e10, 0.0, 1.05], 'Convective zone size - 4.5G & var. rotational velocity', 'cz_var_vel_4_5g');
end
function plot_cz_size_5_0G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_cz_size_plots(gauss_fields(5,:), rotational_vels(2:6,:), true, 0.1, [1.0e2, 1.0e10, 0.0, 1.05], 'Convective zone size - 5.0G & var. rotational velocity', 'cz_var_vel_5_0g');
end
function plot_cz_size_5_5G_var_vel()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_cz_size_plots(gauss_fields(6,:), rotational_vels(2:6,:), true, 0.1, [1.0e2, 1.0e10, 0.0, 1.05], 'Convective zone size - 5.5G & var. rotational velocity', 'cz_var_vel_5_5g');
end

function plot_cz_size_028vc_var_g()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_cz_size_plots(gauss_fields(1:6,:), rotational_vels(5,:), false, 0.1, [1.0e2, 1.0e10, 0.0, 1.05], 'Convective zone radius - vcrit=0.028 & var. magnetic field', 'cz_vc_028_var_g');
end

function plot_cz_size_028vc_var_g_z1()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_cz_size_plots(gauss_fields(1:6,:), rotational_vels(5,:), false, 0.01, [1.0e7, 1.0e10, 0.25, 0.30], 'Convective zone radius - vcrit=0.028 & var. magnetic field','cz_vc_028_var_g_z1');
end

function plot_cz_size_0G_var_vel_z1()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_cz_size_plots(gauss_fields(1,:), rotational_vels(2:6,:), true, 0.01, [1.0e7, 1.0e10, 0.25, 0.30], 'Convective zone radius - 0G & var. rotational velocity', 'cz_var_vel_0_0g_z1');
end



function plot_cz_size_3_5G_var_vel_z1()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_cz_size_plots(gauss_fields(2,:), rotational_vels(2:6,:), true, 0.01, [1.0e7, 1.0e10, 0.25, 0.28], 'Convective zone size - 3.5G & var. rotational velocity', 'cz_var_vel_3_5g_z1');
end


function plot_m_dot_028vc_var_g()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_m_dot_plots(gauss_fields(1:6,:), rotational_vels(5,:), false, 0.5, [1.0e2, 1.0e10, -14.0, -10.0], 'Mass loss - vcrit=0.028 & var. magnetic field', 'mdot_vc_028_var_g');
end

function plot_m_dot_028vc_var_g_z1()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_m_dot_plots(gauss_fields(1:6,:), rotational_vels(5,:), false, 0.05, [3.0e8, 1.0e10, -13.6, -13.3], 'Mass loss - vcrit=0.028 & var. magnetic field', 'mdot_vc_028_var_g_z1');
end



function plot_age_vs_mb_activation_3_5G()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_mb_activation(gauss_fields(2,:), rotational_vels(2:6,:), true, 'Magnetic braking activation & core nature - 3.5G & var. rotational velocity', 'mb_act_var_vel_3_5g');
end

function plot_age_vs_mb_activation_4_0G()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_mb_activation(gauss_fields(3,:), rotational_vels(2:6,:), true, 'Magnetic braking activation & core nature - 4.0G & var. rotational velocity', 'mb_act_var_vel_4_0g');
end

function plot_age_vs_mb_activation_4_5G()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_mb_activation(gauss_fields(4,:), rotational_vels(2:6,:), true, 'Magnetic braking activation & core nature - 4.5G & var. rotational velocity', 'mb_act_var_vel_4_5g');
end

function plot_age_vs_mb_activation_5_0G()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_mb_activation(gauss_fields(5,:), rotational_vels(2:6,:), true, 'Magnetic braking activation & core nature - 5.0G & var. rotational velocity', 'mb_act_var_vel_5_0g');
end

function plot_age_vs_mb_activation_5_5G()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_mb_activation(gauss_fields(6,:), rotational_vels(2:6,:), true, 'Magnetic braking activation & core nature - 5.5G & var. rotational velocity', 'mb_act_var_vel_5_5g');
end

function plot_age_vs_mb_activation_0084vc()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_mb_activation(gauss_fields(2:6,:), rotational_vels(2,:), false, 'Magnetic braking activation & core nature - vcrit=0.0084 & var. magnetic field', 'mb_act_vc_0084_var_g');
end

function plot_age_vs_mb_activation_028vc()
  global gauss_fields;
  global rotational_vels;
  
  age_vs_mb_activation(gauss_fields(2:6,:), rotational_vels(5,:), false, 'Magnetic braking activation & core nature - vcrit=0.028 & var. magnetic field', 'mb_act_vc_028_var_g');
end


function paper1()
  plot_0G_var_vel();
  plot_0G_var_vel_z1();
  plot_vel_rot_0G_var_vel();
  %plot_hr_0336vc_var_g();
  plot_hr_0336vc_var_g_z1();
  plot_hr_0G_var_vel_z1();
  plot_4_0G_var_vel();
  plot_vel_rot_4G_var_vel();
  plot_vel_rot_4G_var_vel_z1();
  plot_cz_size_028vc_var_g();
  plot_cz_size_028vc_var_g_z1();
  plot_m_dot_028vc_var_g();
  plot_m_dot_028vc_var_g_z1();
  plot_age_vs_mb_activation_028vc();
  
  %grid Li var_vel
  plot_0G_var_vel();
  plot_3_5G_var_vel();  
  plot_4_0G_var_vel(); 
  plot_4_5G_var_vel(); 
  plot_5_0G_var_vel(); 
  plot_5_5G_var_vel(); 
 
  %grid Li var_b
  plot_0084vc_var_g();
  plot_014vc_var_g();
  plot_0196vc_var_g();
  plot_028vc_var_g();
  plot_0336vc_var_g();
  
  %grid rot vel
  plot_vel_rot_0G_var_vel();
  plot_vel_rot_3_5G_var_vel();
  plot_vel_rot_4G_var_vel();
  plot_vel_rot_4_5G_var_vel();
  plot_vel_rot_5G_var_vel();
  plot_vel_rot_5_5G_var_vel();
  
  %grid mag breaking
  plot_age_vs_mb_activation_3_5G();
  plot_age_vs_mb_activation_4_0G();
  plot_age_vs_mb_activation_4_5G();
  plot_age_vs_mb_activation_5_0G();
  plot_age_vs_mb_activation_5_5G();
  
end


function main()
  %plot_age_vs_mb_activation_3_5G();
  %plot_age_vs_mb_activation_4_0G();
  %plot_age_vs_mb_activation_4_5G();
  %plot_age_vs_mb_activation_5_0G();
  %plot_age_vs_mb_activation_5_5G();
  %plot_age_vs_mb_activation_0084vc();
  %plot_age_vs_mb_activation_028vc();
  
  
  %plot_0G_var_vel();
  %plot_0G_var_vel_z1();
  %plot_3_5G_var_vel();  
  %plot_4_0G_var_vel(); 
  %plot_4_5G_var_vel(); 
  %plot_5_0G_var_vel(); 
  %plot_5_5G_var_vel(); 
  
  %plot_0084vc_var_g();
  %plot_014vc_var_g();
  %plot_0196vc_var_g();
  %plot_028vc_var_g();
  %plot_0336vc_var_g();
  
  
  %plot_vel_rot_0G_var_vel();
  %plot_vel_rot_3_5G_var_vel();
  %plot_vel_rot_4G_var_vel();
  %plot_vel_rot_4G_var_vel_z1();
  %plot_vel_rot_4_5G_var_vel();
  %plot_vel_rot_5G_var_vel();
  %plot_vel_rot_5_5G_var_vel();  
  
  %plot_cz_size_0G_var_vel();
  %plot_cz_size_3_5G_var_vel();  
  %plot_cz_size_3_5G_var_vel_z_1();
  %plot_cz_size_4_0G_var_vel();  
  %plot_cz_size_4_5G_var_vel();  
  %plot_cz_size_5_0G_var_vel();  
  %plot_cz_size_5_5G_var_vel();  
  %plot_cz_size_028vc_var_g();
  %plot_cz_size_028vc_var_g_z1();
  %plot_cz_size_0G_var_vel_z1();
  
  %plot_m_dot_028vc_var_g();
  %plot_m_dot_028vc_var_g_z1();
    
  
  %plot_hr_3_5G_var_vel();
  %plot_hr_5_0G_var_vel();
  %plot_hr_0336vc_var_g();
  %plot_hr_0336vc_var_g_z1();
  %plot_hr_0G_var_vel();
  %plot_hr_0G_var_vel_z1();
  %plot_hr_3_5G_var_vel_z_1();
  %filename  = "/home/rcaballeron/MESA/workspace/sun-jupiter-system/Docs/runs/run_paper/4g_12kms/1M_photosphere_history.data";
  %calculate_ZAMS(filename);
  
  paper1();
end


main();