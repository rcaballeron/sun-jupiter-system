format longEng;
% grid on by default. Each plot decides if any axis should be hidden
%grid on;
clear global;

global header_lines       = 6;
global table_header_lines = 1;

%cols position in output MESA files
global star_age_col       = 3;
global star_mass_col      = 5;
global log_abs_m_dot_col  = 7;
global mass_conv_core_col = 8;
global log_LH_col         = 17;
global log_Teff_col       = 33;
global log_L_col          = 34;
global log_R_col          = 35;
global log_g_col          = 36;
global surf_avg_omega_col = 39;
global surf_avg_v_rot_col = 42;
#global surf_avg_omega_col = 160;
#global surf_avg_v_rot_col = 163;

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
global sz_top_omega_col   = 87;
global sz_bot_omega_col   = 88;
global core_top_radius_col= 90;
global core_bot_radius_col= 91;
global core_top_omega_col = 98;
global alpha_mlt_col      = 100;
global tau_c_col          = 103;
global rossby_col         = 104;
global rossby_norm_col    = 105;
global b_equi_col         = 106;
global b_star_col         = 107;
global f_min_col          = 108;
global f_max_col          = 109;
global f_star_col         = 110;
global bf_min_col         = 111;
global bf_max_col         = 112;
global bf_star_col        = 113;
global alfven_r_col       = 114;

global num_mix_regions    = 0;
global mix_type_ini_col   = 17;
global mix_relr_top_ini_col = 78;

#Clusters section
global clusters_object_col1 = 1;
global clusters_object_col2 = 2;
global clusters_teff_col = 3;
global clusters_log_g_col = 5;
global clusters_e_log_g_col = 6;
global clusters_fe_h_col = 7;
global clusters_li_col = 9;
global clusters_e_li_col = 10;
global clusters_age_col = 11;
global pleiades_table_filename = 'Pleiades_J_A+A_613_A63_tableb1';
global ori_table_filename = '25_Ori_J_A+A_659_A85_tableb1';
global gamma_vel_a_table_filename = 'Gamma_Vel_A_J_A+A_659_A85_tableb2';
global gamma_vel_b_table_filename = 'Gamma_Vel_B_J_A+A_659_A85_tableb3';
global ngc_2451_b_table_filename = 'NGC_2451_B_J_A+A_659_A85_tableb4';
global ngc_2547_table_filename = 'NGC_2547_J_A+A_659_A85_tableb5';
global ngc_2516_table_filename = 'NGC_2516_J_A+A_659_A85_tableb6';
global ges_dr5_all_oc_table_filename = 'ges_dr5_oc_li1_feh_teff_logg_age2.csv';
global br_21 ='Br21';
global br_22 ='Br22';
#global br_25 ='Br25';
global br_30 ='Br30';
global br_31 ='Br31';
global br_32 ='Br32';
global br_36 ='Br36';
global br_39 ='Br39';
global br_44 ='Br44';
global br_73 ='Br73';
global br_75 ='Br75';
global br_81 ='Br81';
global blanco_1 = 'Blanco1';
global assc_50 = 'Assc50';
global ori_25 = '25_Ori';
global cha_i = 'Cha_I';
global col_197 = 'Col197';
global cz_24 = 'Cz24';
global cz_30 = 'Cz30';
global es_092_05 = 'ES092_05';
global haf_10 = 'Haf10';
global ic_2391 = 'IC2391';
global ic_2602 = 'IC2602';
global ic_4665 = 'IC4665';
#global loden165 = 'Loden165';
global m_67 = 'M67';
global ngc_2141 = 'NGC2141';
global ngc_2158 = 'NGC2158';
global ngc_2232 = 'NGC2232';
global ngc_2243 = 'NGC2243';
global ngc_2355 = 'NGC2355';
global ngc_2244 = 'NGC2244';
global ngc_2264 = 'NGC2264';
global ngc_2420 = 'NGC2420';
global ngc_2425 = 'NGC2425';
global ngc_2451 = 'NGC2451';
global ngc_2516 = 'NGC2516';
global ngc_3532 = 'NGC3532';
global ngc_4815 = 'NGC4815';
global ngc_6005 = 'NGC6005';
global ngc_6067 = 'NGC6067';
global ngc_6253 = 'NGC6253';
global ngc_6259 = 'NGC6259';
global ngc_6281 = 'NGC6281';
global ngc_6405 = 'NGC6405';
global ngc_6530 = 'NGC6530';
global ngc_6633 = 'NGC6633';
global ngc_6649 = 'NGC6649';
global ngc_6705 = 'NGC6705';
global ngc_6709 = 'NGC6709';
global ngc_6802 = 'NGC6802';
global pismis_15 = 'Pismis15';
global pismis_18 = 'Pismis18';
#global rho_oph = 'Rho_Oph';
global rup_134 = 'Rup134';
global trumpler_14 = 'Trumpler14';
global trumpler_20 = 'Trumpler20';
global trumpler_23 = 'Trumpler23';
global trumpler_5 = 'Trumpler5';
global gamma2_vel = 'gamma2_Vel';
global lam_ori = 'lam_Ori';

global teff_low_limit = 5700;
global teff_top_limit = 5800;



global data_parent_folder = '/home/rcaballeron/MESA/workspace/sun-jupiter-system/Docs/runs/run_thesis';
global tables_parent_folder = '/home/rcaballeron/MESA/workspace/sun-jupiter-system/project/tables';
global filename = '1M_photosphere_history.data';

%global filename = 'annexB/1M_photosphere_history.data';
%global filename = 'annex2B/1M_photosphere_history.data';
global gauss_fields = ['0g'; '2g'; '2.5g'; '3g'; '3.3g'; '3.5g'; '4g'; '4.3g'; '4.5g'; '5g'; '5.5g'; 'Xg'];
%global gauss_fields = ['0g'; '3g'; '3.5g'; '4g'; '4.5g'; '5g'];
global rotational_vels = ['0crit';'0084crit'; '014crit'; '0196crit'; '028crit';
  '0336crit';'029crit';'030crit';'031crit';'0312crit';'0314crit';'032crit';'0336crit_alpha';
  '11crit';'105crit';'1075crit';'1125crit';'1025crit';
  '115crit';'10crit';'0975crit';'095crit';'0925crit';'125crit';
  '135crit';'1175crit';'12crit';'1225crit';'1275crit';
  '13crit';'1325crit';'1375crit';'14crit';'1425crit';
  '145crit';'1475crit';'15crit';'1525crit';'155crit'];
%dl -> disk locking
global dl_rotational_vels = ['0crit_dl';'0084crit_dl'; '014crit_dl'; '0196crit_dl';
  '028crit_dl'; '0336crit_dl';'029crit_dl';'030crit_dl';'031crit_dl';
  '0312crit_dl';'0314crit_dl';'032crit_dl';'9_090256e-6_dl'];
global colors = ['k'; 'r'; 'g'; 'b'; 'y'; 'm'; 'c'];

%Array indexes
global idx_0_0G  = 1;
global idx_2_0G  = 2;
global idx_2_5G  = 3;
global idx_3_0G  = 4;
global idx_3_3G  = 5;
global idx_3_5G  = 6;
global idx_4_0G  = 7;
global idx_4_3G  = 8;
global idx_4_5G  = 9;
global idx_5_0G  = 10;
global idx_5_5G  = 11;
global idx_X_G  = 12;

%The following are with mlt=1.82
global idx_0crit    = 1;
global idx_0084crit = 2;
global idx_014crit  = 3;
global idx_0196crit = 4;
global idx_028crit  = 5;
global idx_0336crit = 6;
%The following are with mlt=1.7
global idx_029crit  = 7;
global idx_030crit  = 8;
global idx_031crit  = 9;
global idx_0312crit  = 10;
global idx_0314crit  = 11;
global idx_032crit  = 12;
global idx_0336crit_alpha = 13;
global idx_11crit  = 14;
global idx_105crit  = 15;
global idx_1075crit  = 16;
global idx_1125crit  = 17;
global idx_1025crit  = 18;
global idx_115crit  = 19;
global idx_10crit  = 20;
global idx_0975crit  = 21;
global idx_095crit  = 22;
global idx_0925crit  = 23;
global idx_125crit  = 24;
global idx_135crit  = 25;
global idx_1175crit  = 26;
global idx_12crit  = 27;
global idx_1225crit  = 28;
global idx_1275crit  = 29;
global idx_13crit  = 30;
global idx_1325crit  = 31;
global idx_1375crit  = 32;
global idx_14crit  = 33;
global idx_1425crit  = 34;
global idx_145crit  = 35;
global idx_1475crit  = 36;
global idx_15crit  = 37;
global idx_1525crit  = 38;
global idx_155crit  = 39;


%The following with mlt=var and disk locking
global idx_9_090256e_6_dl  = 13;



%Sun constants
global sun_gauss_field = 1.0; %G
global sun_age = 4.57e9; %years
global sun_A_Li7 = 1.1;
global sun_A_eLi7 = 0.1;
global sun_vel_rot = 2.0; %km/s
global sun_T_eff = 5772; %K
global sun_L = 3.828e26; %W
global sun_omega = 2.903e-6 %rad/s
global sun_log_g = 4.4374
global sun_radius = 6.957e8

%Pleiades constants
%global pleiades_age = 1.0e8; %years
global pleiades_A_Li7 = 2.95;


%ZAMS
%global delta_h1_center = 0.0015;
global delta_h1_center = 0.0001; %When max_time fixed to 1e6 years
global ratio_log_LH_vs_log_L = 0.99;

%Graphic size in inches
global line_width = 3;
global W = 18;
global H = 14;
global tick_font_size = 32;
global title_font_size = 32;
global axis_font_size = 32;
global legend_font_size = 27;

%For the axis
global gigaYear = 1000000000;





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
  filename

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
  A = [A zeros(size(A, 1), 2)];

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
  %delta_h1_center
  %A
  %Find those row in which diff h1 > delta_h1_center
  B = find(A(:,5) > delta_h1_center);
  %Select the first one and get rows from A till that row
  C = A(1:B(1),:);
  % Now find those row in which LH div L > ratio_log_LH_vs_log_L
  D = find(C(:,6) > ratio_log_LH_vs_log_L);

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

function [plot_path] = plot_clusters(A, color, width, ytick, axis_limits, aidx, subfolder)
  global tables_parent_folder;
  global pleiades_table_filename;
  global ori_table_filename;
  global gamma_vel_a_table_filename;
  global gamma_vel_b_table_filename;
  global ngc_2451_b_table_filename;
  global ngc_2547_table_filename;
  global ngc_2516_table_filename;
  global ges_dr5_all_oc_table_filename;
  global br_21;
  global br_22;
  #global br_25;
  global br_30;
  global br_31;
  global br_32;
  global br_36;
  global br_39;
  global br_44;
  global br_73;
  global br_75;
  global br_81;
  global blanco_1;
  global assc_50;
  global ori_25;
  global cha_i;
  global col_197;
  global cz_24;
  global cz_30;
  global es_092_05;
  global haf_10;
  global ic_2391;
  global ic_2602;
  global ic_4665;
  #global loden165;
  global m_67;
  global ngc_2141;
  global ngc_2158;
  global ngc_2232;
  global ngc_2243;
  global ngc_2355;
  global ngc_2244;
  global ngc_2264;
  global ngc_2420;
  global ngc_2425;
  global ngc_2451;
  global ngc_2516;
  global ngc_3532;
  global ngc_4815;
  global ngc_6005;
  global ngc_6067;
  global ngc_6253;
  global ngc_6259;
  global ngc_6281;
  global ngc_6405;
  global ngc_6530;
  global ngc_6633;
  global ngc_6649;
  global ngc_6705;
  global ngc_6709;
  global ngc_6802;
  global pismis_15;
  global pismis_18;
  #global rho_oph;
  global rup_134;
  global trumpler_14;
  global trumpler_20;
  global trumpler_23;
  global trumpler_5;
  global gamma2_vel;
  global lam_ori;




 %{
  full_path = strcat(tables_parent_folder, '/', pleiades_table_filename);
  plot_cluster(A, 'k', '+', width, ytick, axis_limits, full_path);

  full_path = strcat(tables_parent_folder, '/', ori_table_filename);
  plot_cluster(A, 'r', 'o', width, ytick, axis_limits, full_path);

  full_path = strcat(tables_parent_folder, '/', gamma_vel_a_table_filename);
  plot_cluster(A, 'g', '*', width, ytick, axis_limits, full_path);

  full_path = strcat(tables_parent_folder, '/', gamma_vel_b_table_filename);
  plot_cluster(A, 'b', 'x', width, ytick, axis_limits, full_path);

  full_path = strcat(tables_parent_folder, '/', ngc_2451_b_table_filename);
  plot_cluster(A, 'y', 's', width, ytick, axis_limits, full_path);

  full_path = strcat(tables_parent_folder, '/', ngc_2547_table_filename);
  plot_cluster(A, 'm', 'd', width, ytick, axis_limits, full_path);

  full_path = strcat(tables_parent_folder, '/', ngc_2516_table_filename);
  plot_cluster(A, 'c', '^', width, ytick, axis_limits, full_path);

  %}

  full_path = strcat(tables_parent_folder, '/', br_21);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', br_22);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  #full_path = strcat(tables_parent_folder, '/', br_25);
  #plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', br_30);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', br_31);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', br_32);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', br_36);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', br_39);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', br_44);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', br_73);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', br_75);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', br_81);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', cha_i);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', col_197);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', cz_24);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', cz_30);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', es_092_05);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', haf_10);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ic_2391);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ic_2602);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ic_4665);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  #full_path = strcat(tables_parent_folder, '/', loden165);
  #plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', m_67);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_2141);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_2158);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_2232);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_2243);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_2244);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_2264);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_2355);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_2420);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_2425);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_2451);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_2516);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_3532);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_4815);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_6005);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_6067);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_6253);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_6259);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_6281);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_6405);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_6530);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_6633);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_6649);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_6705);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_6709);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', ngc_6802);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', pismis_15);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', pismis_18);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  #full_path = strcat(tables_parent_folder, '/', rho_oph);
  #plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', rup_134);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', trumpler_14);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', trumpler_20);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', trumpler_23);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', trumpler_5);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', gamma2_vel);
  plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

  full_path = strcat(tables_parent_folder, '/', lam_ori);
  plot_path = plot_cluster(A, color, '+', width, ytick, axis_limits, full_path, aidx, subfolder);

%{
  full_path = strcat(tables_parent_folder, '/', blanco_1);
  plot_cluster(A, 'r', 'o', width, ytick, axis_limits, full_path);

  full_path = strcat(tables_parent_folder, '/', assc_50);
  plot_cluster(A, 'g', '*', width, ytick, axis_limits, full_path);

  full_path = strcat(tables_parent_folder, '/', ori_25);
  plot_cluster(A, 'b', 'x', width, ytick, axis_limits, full_path);
%}

end

function [root_subfolder] = plot_cluster(A, color, marker, width, ytick, axis_limits, full_path, aidx, subfolder)
  global tables_parent_folder
  global tick_font_size;
  global table_header_lines;
  global clusters_object_col1;
  global clusters_object_col2;
  global clusters_teff_col;
  global clusters_li_col;
  global clusters_e_li_col;
  global clusters_age_col;
  global clusters_log_g_col;
  global clusters_e_log_g_col;
  global clusters_fe_h_col
  global gigaYear;
  delta_log_g = 0.05 %dex
  delta_teff = 50 %K
  delta_fe_h = 0.05 %dex

  %fmt = get_parsing_fmt([clusters_teff_col, clusters_li_col, clusters_age_col, clusters_err_up_age_col, clusters_err_down_age_col]);
  %fmt = get_parsing_fmt([clusters_teff_col, clusters_log_g_col, clusters_e_log_g_col, clusters_li_col, clusters_age_col]);
  fmt = get_parsing_fmt([clusters_object_col1, clusters_object_col2, clusters_teff_col, clusters_log_g_col, clusters_fe_h_col, clusters_li_col, clusters_e_li_col, clusters_age_col]);

  %C = read_matrix_from_file(full_path, fmt, table_header_lines, 5);
  C = read_matrix_from_file(full_path, fmt, table_header_lines, 8);


  %Get cluster age (in yrs) and limits
  cluster_age = C(1, 8) * gigaYear;
  %cluster_top_age = C(1, 4) + cluster_age;
  %cluster_low_age = C(1, 5) + cluster_age;
  % Overall, the uncertainty on the determination of log t ranges from 0.15 to 0.25
  % for young clusters (<1Gyr) and from 0.1 to 0.2 for old clusters.
  if (cluster_age >= gigaYear)
    age_top_limit = cluster_age + ((cluster_age*0.1)/100);
    age_low_limit = cluster_age - ((cluster_age*0.1)/100);
  else
    age_top_limit = cluster_age + ((cluster_age*0.15)/100);
    age_low_limit = cluster_age - ((cluster_age*0.15)/100);
  endif

  %Find rows with age older than cluster_low_age and get: Teff, logg
  D = find(A(:,1) >= age_low_limit);
  %teff is given in log10 and is present in A(2)
  teff_low_limit = power(10, A(D(1), 2)) - delta_teff;
  %logg is present in A(3)
  log_g_low_limit = A(D(1), 3) - delta_log_g;
  %FeH compared to the Sun one
  fe_h_low_limit = 0.0 - delta_fe_h;

  %Find rows with age older than cluster_top_age and get Teff, logg
  D = find(A(:,1) >= age_top_limit);
  teff_top_limit = power(10, A(D(1), 2)) + delta_teff;
  log_g_top_limit = A(D(1), 3) + delta_log_g;
  fe_h_top_limit = 0.0 + delta_fe_h;
  %{
  cluster_top_age
  cluster_low_age
  teff_top_limit
  teff_low_limit
  log_g_low_limit
  log_g_top_limit
%}
  %First filter by Teff
  %Filter out starts with Teff below and above the limits
  B1 = C(:,3) < teff_low_limit;
  B2 = C(:,3) > teff_top_limit;

  % Mark rows which fulfill either B1 or B2 and remove them
  B = B1 | B2;
  C(B,:) = [];
  %C

  %Second filter by logg
  B1 = C(:,4) < log_g_low_limit;
  B2 = C(:,4) > log_g_top_limit;
  B = B1 | B2;
  C(B,:) = [];
  %C

  %Third filter by metallicity
  B1 = C(:,5) < fe_h_low_limit;
  B2 = C(:,5) > fe_h_top_limit;
  B = B1 | B2;
  C(B,:) = [];


  %Create subfolder for each picture and inside it for each data file
  folder = substr(full_path, 1, rindex(full_path, "/"));
  name = substr(full_path, rindex(full_path, "/") + 1);
  root_subfolder = strcat(folder, aidx, "/");
  new_subfolder = strcat(root_subfolder, subfolder);
  mkdir(new_subfolder);
  filename = strcat(new_subfolder, "/", name, ".sub");
  fid = fopen (filename, "w");



  filter_setup = strcat("%Filter setup: $\\teff$ low=", num2str(teff_low_limit), ", $\\teff$ top=", num2str(teff_top_limit),
    ", $\\gsurf$ low=", num2str(log_g_low_limit), ", $\\gsurf$ top=", num2str(log_g_top_limit),
    ", $\\feh$ low=", num2str(fe_h_low_limit), ", $\\feh$ top=", num2str(fe_h_top_limit),
    ", Age low=", num2str(age_low_limit/1000000000), ", Age top=", num2str(age_top_limit/1000000000),
    "\n");

  fputs (fid, filter_setup);
  fputs (fid, "ObjectId Teff(K) logg(dex) FeH(dex) ALi(dex) eALi(dex) Age(Gyr)\n");

  fclose (fid);
  dlmwrite (filename, C, "delimiter", " & ", "newline", "\\\\\\\\ \n", "-append");


  %Plot values
  %plot(A(:,3), A(:,2), marker, 'markersize', 8, 'color', [0.5,0.1,0.8], 'markerfacecolor', [0.5,0.1,0.8]);
  %plot(C(:,6)*1000000000, C(:,4), marker, 'markersize', 8, 'color', color, 'markerfacecolor', color);
  %h = errorbar(C(:,6)*1000000000, C(:,4),C(:,5),"~d");
  h = errorbar(C(:,8)*1000000000, C(:,6),C(:,7),"~d");
  set(h, 'markersize', 5, 'color', color, 'markerfacecolor', color, 'linewidth', 1.5);
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

  A_Li7(:,1) = fA_Li7(A(:,4),A(:,5));

  result = A_Li7;
end


function plot_A_Li7(A, B, color, width, ytick, axis_limits)
  global tick_font_size;
  global gigaYear;

  %Plot values
  plot(A(:,1), B(:,1), color, 'linewidth', width);

  %Axis scales
  grid on;
  set(gca, 'XScale', 'log');
  set(gca,'xminorgrid','off');

  %Axis limits
  axis(axis_limits);

  %Axis ticks
  %set(gca,'YTick',1.0:ytick:4.5);

  xticks = get (gca, "xtick");
  %xlabels = arrayfun (@(x) sprintf ("%.2e", x), xticks, "uniformoutput", false);
  %User Ga for the age instead of scientific notation
  %xlabels = arrayfun (@(x) sprintf ("%2.3f", x/gigaYear), xticks, "uniformoutput", false);
  xlabels = arrayfun (@(x) sprintf ("%g", x/gigaYear), xticks, "uniformoutput", false);
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
  grid on;
  set(gca,'xminorgrid','off');
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

  %
  %errorbar (x, y1, err, err, err, err, "~>"

end



function plot_vel_rot(A, color, width, ytick, axis_limits)
  global tick_font_size;
  global gigaYear;

  %Plot values
  plot(A(:,1), A(:,2), color, 'linewidth', width);
  %Don't show sup limit of CZ
  %plot(A(:,1), A(:,3), color, 'linewidth', width, 'linestyle', '--');
  plot(A(:,1), A(:,4), color, 'linewidth', width, 'linestyle', ':');

  %Axis scales
  set(gca, 'XScale', 'log');
  grid on;
  set(gca,'xminorgrid','off');


  set(gca,'YTick',axis_limits(3):ytick:axis_limits(4));
  %set(gca,'YTick',0:2.5:25);
  %set(gca,'XTick',axis_limits(1):10:axis_limits(2));
  %set(gca,'XTick',[1.0e5,1.0e6,1.0e10]);
  %set(gca,'XTickLabel',[])


  %Axis ticks
  axis(axis_limits);
  %axis([axis_limits(1), axis_limits(2), 0, 20]);
  xticks = get (gca, "xtick");
  %xlabels = arrayfun (@(x) sprintf ("%.1e", x), xticks, "uniformoutput", false);
  %User Ga for the age instead of scientific notation
  xlabels = arrayfun (@(x) sprintf ("%g", x/gigaYear), xticks, "uniformoutput", false);
  set (gca, "xticklabel", xlabels) ;
  %set(gca, "xminortick","on", "xminorgrid","off")
  yticks = get (gca, "ytick");
  ylabels = arrayfun (@(x) sprintf ("%2.1f", x), yticks, "uniformoutput", false);
  set (gca, "yticklabel", ylabels);
  set(gca, 'fontsize', tick_font_size);

end




function plot_omega(A, color, width, ytick, axis_limits)
  global tick_font_size;
  global sun_omega;
  global gigaYear;

  %Plot values
  plot(A(:,1), A(:,2) ./ sun_omega, color, 'linewidth', width);
  plot(A(:,1), A(:,3) ./ sun_omega, 'linewidth', width, 'linestyle', ':');

  %Axis scales
  set(gca, 'XScale', 'log');
  set(gca,'YScale', 'log');
  grid on;
  set(gca,'xminorgrid','off');


  %Axis ticks
  axis(axis_limits);
  xticks = get (gca, "xtick");
  %xlabels = arrayfun (@(x) sprintf ("%.2e", x), xticks, "uniformoutput", false);
  %User Ga for the age instead of scientific notation
  xlabels = arrayfun (@(x) sprintf ("%g", x/gigaYear), xticks, "uniformoutput", false);
  set (gca, "xticklabel", xlabels) ;
  yticks = get (gca, "ytick");

  ylabels = arrayfun (@(x) sprintf ("%2.2f", x), yticks, "uniformoutput", false);
  set (gca, "yticklabel", ylabels);
  set(gca, 'fontsize', tick_font_size);
end

function plot_teff(A, color, width, ytick, axis_limits)
  global tick_font_size;
  global sun_omega;
  global gigaYear;


  %Plot values
  plot(A(:,1), A(:,2), color, 'linewidth', width);
  plot(A(:,1), A(:,3), color, 'linewidth', width, 'linestyle', ':');

  %Axis scales
  set(gca, 'XScale', 'log');
  grid on;
  set(gca,'xminorgrid','off');

  %Axis ticks
  axis(axis_limits);

  xticks = get (gca, "xtick");
  %xlabels = arrayfun (@(x) sprintf ("%.2e", x), xticks, "uniformoutput", false);
  %User Ga for the age instead of scientific notation
  xlabels = arrayfun (@(x) sprintf ("%g", x/gigaYear), xticks, "uniformoutput", false);
  set (gca, "xticklabel", xlabels) ;

  yticks = get (gca, "ytick");
  ylabels = arrayfun (@(x) sprintf ("%2.1f", x), yticks, "uniformoutput", false);
  set (gca, "yticklabel", ylabels);

  set(gca, 'fontsize', tick_font_size);


end


function plot_radius(A, color, width, ytick, axis_limits)
  global tick_font_size;
  global sun_omega;
  global gigaYear;

  %Plot values
  plot(A(:,1), A(:,2), color, 'linewidth', width);

  %Axis scales
  set(gca, 'XScale', 'log');
  grid on;
  set(gca,'xminorgrid','off');

  %Axis ticks
  axis(axis_limits);

  xticks = get (gca, "xtick");
  %xlabels = arrayfun (@(x) sprintf ("%.2e", x), xticks, "uniformoutput", false);
  %User Ga for the age instead of scientific notation
  xlabels = arrayfun (@(x) sprintf ("%g", x/gigaYear), xticks, "uniformoutput", false);
  set (gca, "xticklabel", xlabels) ;

  yticks = get (gca, "ytick");
  ylabels = arrayfun (@(x) sprintf ("%2.1f", x), yticks, "uniformoutput", false);
  set (gca, "yticklabel", ylabels);

  set(gca, 'fontsize', tick_font_size);

end



function plot_size_cz(A, color, width, ytick, axis_limits)
  plot_size_cz2(A, color, width, ytick, axis_limits, 0, '-');
end

function plot_size_cz2(A, color, width, ytick, axis_limits, offset, linestyle)
  global tick_font_size;
  global gigaYear;

  %Plot values
  %plot(A(:,1), A(:,3) .- A(:,4), color, 'linewidth', width);
  plot(A(:,1), A(:,3) .- A(:,4) .+ offset , color, 'linewidth', width, 'linestyle', linestyle);
  %plot(A(:,1), A(:,5) .- A(:,6) .+ offset , color, 'linewidth', width, 'linestyle', '--');
  %plot(A(:,1), 10 .** A(:,2), color, 'linewidth', width, 'linestyle', '--');

  %Axis scales
  set(gca, 'XScale', 'log');

  %Axis limits
  set(gca,'YTick',0:ytick:1.0);
  grid on;
  set(gca,'xminorgrid','off');

  %Axis ticks
  axis(axis_limits);
  xticks = get (gca, "xtick");
  %xlabels = arrayfun (@(x) sprintf ("%.2e", x), xticks, "uniformoutput", false);
  %User Ga for the age instead of scientific notation
  xlabels = arrayfun (@(x) sprintf ("%g", x/gigaYear), xticks, "uniformoutput", false);
  %Usar la siguiente instrucción si queremos fijar nosotros las marcas a mostrar
  %No está pensado usarlas para las gráficas de zoom
  %xlabels = ['1.00e+02';' ';'1.00e+04';' ';'1.00e+06';' ';'1.00e+08';' ';'1.00e+10'];
  set (gca, "xticklabel", xlabels);
  yticks = get (gca, "ytick");
  ylabels = arrayfun (@(x) sprintf ("%1.2f", x), yticks, "uniformoutput", false);
  set (gca, "yticklabel", ylabels);
  set(gca, 'fontsize', tick_font_size);

end

function plot_m_dot(A, color, width, ytick, axis_limits)
  global tick_font_size;
  global gigaYear;

  %Plot values
  plot(A(:,1), A(:,2), color, 'linewidth', width);

  %Axis scales
  set(gca, 'XScale', 'log');

  %Axis limits
  set(gca,'YTick',axis_limits(3):ytick:axis_limits(4));
  grid on;
  set(gca,'xminorgrid','off');

  %Axis ticks
  axis(axis_limits);
  xticks = get (gca, "xtick");
  %xlabels = arrayfun (@(x) sprintf ("%.2e", x), xticks, "uniformoutput", false);
  xlabels = arrayfun (@(x) sprintf ("%g", x/gigaYear), xticks, "uniformoutput", false);
  %xlabels = ['1.00e+02';' ';'1.00e+04';' ';'1.00e+06';' ';'1.00e+08';' ';'1.00e+10'];
  set (gca, "xticklabel", xlabels) ;
  yticks = get (gca, "ytick");
  ylabels = arrayfun (@(x) sprintf ("%1.2f", x), yticks, "uniformoutput", false);
  set (gca, "yticklabel", ylabels);
  set(gca, 'fontsize', tick_font_size);

end

function plot_reimers_m_dot(A, color, width, ytick, axis_limits)
  global tick_font_size;

  %Plot values
  plot(A(:,1), (power(10,A(:,3)) .* power(10,A(:,4))) ./ A(:,2), color, 'linewidth', width);

  %Axis scales
  set(gca, 'XScale', 'log');

  %Axis limits
  set(gca,'YTick',axis_limits(3):ytick:axis_limits(4));
  grid on;

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

function plot_m_l_r(A, color, width, ytick, axis_limits)
  global tick_font_size;

  %Plot values
  plot(A(:,1), A(:,2), color, 'linewidth', width);
  plot(A(:,1), power(10,A(:,3)), color, 'linewidth', width, 'linestyle', '--');
  plot(A(:,1), power(10,A(:,4)), color, 'linewidth', width, 'linestyle', ':');

  %Axis scales
  set(gca, 'XScale', 'log');

  %Axis limits
  set(gca,'YTick',axis_limits(3):ytick:axis_limits(4));
  grid on;

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

function plot_alpha_mlt(A, color, width, ytick, axis_limits)
  global tick_font_size;
  global gigaYear;

  %Plot values
  plot(A(:,1), A(:,2), color, 'linewidth', width);

  %Axis scales
  set(gca, 'XScale', 'log');

  set(gca,'YTick',axis_limits(3):ytick:axis_limits(4));
  grid on;
  set(gca,'xminorgrid','off');

  %Axis limits
  axis(axis_limits);

  %Axis ticks
  %set(gca,'YTick',1.0:ytick:4.5);

  xticks = get (gca, "xtick");
  %xlabels = arrayfun (@(x) sprintf ("%.2e", x), xticks, "uniformoutput", false);
  %User Ga for the age instead of scientific notation
  xlabels = arrayfun (@(x) sprintf ("%g", x/gigaYear), xticks, "uniformoutput", false);
  set (gca, "xticklabel", xlabels) ;

  yticks = get (gca, "ytick");
  ylabels = arrayfun (@(x) sprintf ("%2.2f", x), yticks, "uniformoutput", false);
  set (gca, "yticklabel", ylabels);

  set(gca, 'fontsize', tick_font_size);

end

function plot_mag_field(A, color, width, show_limits, ytick, axis_limits)
  global tick_font_size;
  global sun_omega;
  global gigaYear;

  %Plot values
  %Omega
  plot(A(:,1), A(:,2) ./ sun_omega, color, 'linewidth', width, 'linestyle', ':');
  if (show_limits)
    %Min bf
    plot(A(:,1), A(:,3), color, 'linewidth', width-2, 'linestyle', '--');
    %plot(A(:,1), A(:,3), color+1, 'linewidth', width, 'linestyle', '--');
    %Max bf
    plot(A(:,1), A(:,4), color, 'linewidth', width-2, 'linestyle', '--');
    %plot(A(:,1), A(:,4), color+2, 'linewidth', width, 'linestyle', '--');
  endif
  %Mean bf
  plot(A(:,1), A(:,5), color, 'linewidth', width);
  %plot(A(:,1), A(:,5), color+3, 'linewidth', width);

  %Axis scales
  set(gca, 'XScale', 'log');
  set(gca, 'YScale', 'log');
  grid on;
  set(gca,'xminorgrid','off');

  %Axis ticks
  axis(axis_limits);
  xticks = get (gca, "xtick");
  %xlabels = arrayfun (@(x) sprintf ("%.1e", x), xticks, "uniformoutput", false);
  %User Ga for the age instead of scientific notation
  xlabels = arrayfun (@(x) sprintf ("%g", x/gigaYear), xticks, "uniformoutput", false);
  set (gca, "xticklabel", xlabels) ;
  yticks = get (gca, "ytick");
  %ylabels = arrayfun (@(x) sprintf ("%2.1f", x), yticks, "uniformoutput", false);
  ylabels = arrayfun (@(x) sprintf ("%2.0f", x), yticks, "uniformoutput", false);
  set (gca, "yticklabel", ylabels);
  set(gca, 'fontsize', tick_font_size);

end


function plot_kipperhahn(A, B, color, width, ytick, axis_limits)
  global tick_font_size;
  global num_mix_regions;
  global colors;

  %Axis scales
  set(gca, 'XScale', 'log');

  %Axis ticks
  axis(axis_limits);
  xticks = get (gca, "xtick");
  xlabels = arrayfun (@(x) sprintf ("%.2e", x), xticks, "uniformoutput", false);
  set (gca, "xticklabel", xlabels) ;
  yticks = get (gca, "ytick");
  ylabels = arrayfun (@(x) sprintf ("%1.1f", x), yticks, "uniformoutput", false);
  set (gca, "yticklabel", ylabels);
  set(gca, 'fontsize', tick_font_size);


  %Plot values
  for i=1:rows(A)
      if (mod (i, 100) == 0)
        i
      endif

    for j=2:columns(A)
      mix_type = A(i,j);

      if (mix_type == 0)
          color = 1;
      else
          color = 2;
      endif


      if (j==2)
        y0 = 0.0;
      else
        y0 = B(i,j-1);
      endif
      y1 = B(i,j);

      line("xdata",[A(i,1),A(i,1)], "ydata",[y0,y1], "linewidth", 1, "color", colors(color,:));
      %drawLine([A(i,1) y0  A(i,1) y1], "linewidth", 1, "color", colors(color,:));
    end
  end
end




function hr_plots(gauss_fields, rotational_vels, is_var_vel, xy_ticks, axis_limits, leg_loc, atitle, afilename, aidx)
  global data_parent_folder;
  global tables_parent_folder;
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
      full_path

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

  %l = legend(labels, "location", "southoutside", "orientation", "horizontal");
  l = legend(labels, "location", leg_loc);
  set (l, "fontsize", legend_font_size);
  %legend boxoff
  xlabel('log Teff', 'fontsize', axis_font_size);
  ylabel('log (L/L_{sun})', 'fontsize', axis_font_size);
  %ylabel('{\bf \omega} = a {\bf V}');
  title(atitle, 'fontsize', title_font_size);

  hold('off');

  plot_path = strcat(tables_parent_folder, '/', aidx);
  save_figure(f, strcat(plot_path, '/', afilename, aidx));
end



function age_vs_li_plots(gauss_fields, rotational_vels, is_var_vel, ytick, axis_limits, leg_loc, atitle, afilename, aidx)
  global data_parent_folder;
  global filename;
  global star_age_col;
  global surface_h1_col;
  global surface_li_col;
  global log_Teff_col;
  global log_g_col;
  global header_lines;
  global sun_age;
  global sun_A_Li7;
  global sun_A_eLi7;
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

  %First plot Li
  for i=1:rows(gauss_fields)
    for j=1:rows(rotational_vels)
      sub_folder = strcat(gauss_fields(i,:), '_', rotational_vels(j,:));
      full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);

      fmt = get_parsing_fmt([star_age_col, log_Teff_col, log_g_col, surface_h1_col, surface_li_col]);

      A = read_matrix_from_file(full_path, fmt, header_lines, 5);

      B = calculate_A_Li7(A);

      plot_A_Li7(A, B, colors(i*j,:), line_width, ytick, axis_limits);

      % Plot ZAMS reference
      %zams = calculate_ZAMS(full_path);
      %line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 2, "linestyle", "--", "color", colors(i*j,:))

      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['A(Li)-', strtrim(rotational_vels(j,:))]};
        %labels = {labels{:}, ['ZAMS-', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['A(Li)-', strtrim(gauss_fields(i,:))]};
        %labels = {labels{:}, ['ZAMS-', strtrim(gauss_fields(i,:))]};
      endif

      %plot_clusters(A, line_width, ytick, axis_limits); <-- Doing this the legend is not shown properly

    end
  end

  %Second plot clusters. This approch, althoug slower, allows to show the legend properly
  for i=1:rows(gauss_fields)
    for j=1:rows(rotational_vels)
      sub_folder = strcat(gauss_fields(i,:), '_', rotational_vels(j,:));
      full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);

      fmt = get_parsing_fmt([star_age_col, log_Teff_col, log_g_col, surface_h1_col, surface_li_col]);

      A = read_matrix_from_file(full_path, fmt, header_lines, 5);
      plot_path = plot_clusters(A, colors(i*j,:), line_width, ytick, axis_limits, aidx, sub_folder);
    end
  end

  %Plot ZAMS, only one of them
  sub_folder = strcat(gauss_fields(1,:), '_', rotational_vels(1,:));
  full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
  zams = calculate_ZAMS(full_path);
  line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 3, "linestyle", "--", "color", "k");

  % Plot sun reference
  h = errorbar(sun_age, sun_A_Li7,sun_A_eLi7,"~d");
  set(h, 'markersize', 5, 'color', [0.5,0.1,0.8], 'markerfacecolor', [0.5,0.1,0.8], 'linewidth', 1.5);
  plot(sun_age, sun_A_Li7, '*', 'markersize', 15, 'color', [0.5,0.1,0.8]);


  %plot(pleiades_age, pleiades_A_Li7, 's', 'markersize', 10, 'color', [0.5,0.1,0.8], 'markerfacecolor', [0.5,0.1,0.8]);

  %l = legend(labels, "location", "southeastoutside");
  l = legend(labels, "location", leg_loc);
  set (l, "fontsize", legend_font_size);
  %legend boxoff
  xlabel('star age (Ga)', 'fontsize', axis_font_size);
  ylabel('A(Li7)', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');

  save_figure(f, strcat(plot_path, '/', afilename, aidx));

end

function age_vs_cz_size_plots(gauss_fields, rotational_vels, is_var_vel, ytick, axis_limits, leg_loc, atitle, afilename, aidx)
  global data_parent_folder;
  global tables_parent_folder;
  global filename;
  global star_age_col;
  global log_R_col;
  global sz_top_radius_col;
  global sz_bot_radius_col;
  global core_top_radius_col;
  global core_bot_radius_col;
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

      fmt = get_parsing_fmt([star_age_col, log_R_col, sz_top_radius_col, sz_bot_radius_col, core_top_radius_col, core_bot_radius_col]);

      A = read_matrix_from_file(full_path, fmt, header_lines, 6);

      plot_size_cz(A, colors(i*j,:), line_width, ytick, axis_limits);

      % Plot ZAMS reference
      %zams = calculate_ZAMS(full_path);
      %line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 2, "linestyle", "--", "color", colors(i*j,:))

      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['CZ_{rad}-', strtrim(rotational_vels(j,:))]};
        %labels = {labels{:}, ['Core_{rad}-', strtrim(rotational_vels(j,:))]};
        %labels = {labels{:}, ['ZAMS-', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['CZ_{rad}-', strtrim(gauss_fields(i,:))]};
        %labels = {labels{:}, ['Core_{rad}-', strtrim(gauss_fields(i,:))]};
        %labels = {labels{:}, ['ZAMS-', strtrim(gauss_fields(i,:))]};
      endif
    end
  end

  %Plot ZAMS, only one of them
  sub_folder = strcat(gauss_fields(1,:), '_', rotational_vels(1,:));
  full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
  zams = calculate_ZAMS(full_path);
  line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 3, "linestyle", "--", "color", "k");


  l = legend(labels, "location", leg_loc);
  set (l, "fontsize", legend_font_size);
  %legend boxoff
  xlabel('star age (Ga)', 'fontsize', axis_font_size);
  ylabel('Size conv. zone (R_{cz}/R_{sun})', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');

  plot_path = strcat(tables_parent_folder, '/', aidx);
  save_figure(f, strcat(plot_path, '/', afilename, aidx));
end

function age_vs_teff_plots(gauss_fields, rotational_vels, is_var_vel, ytick, axis_limits, leg_loc, atitle, afilename, aidx)
  global data_parent_folder;
  global tables_parent_folder;
  global filename;
  global star_age_col;
  global log_Teff_col;
  global log_g_col;
  global header_lines;
  global sun_age;
  global sun_vel_rot;
  global colors;
  global sun_T_eff;
  global sun_log_g;
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

      fmt = get_parsing_fmt([star_age_col, log_Teff_col, log_g_col]);

      A = read_matrix_from_file(full_path, fmt, header_lines, 3);

      axis_limits
      plot_teff(A, colors(i*j,:), line_width, ytick, axis_limits);


      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['log T_{eff}-', strtrim(rotational_vels(j,:))]};
        labels = {labels{:}, ['log g-', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['log T_{eff}-', strtrim(gauss_fields(i,:))]};
        labels = {labels{:}, ['log g-', strtrim(gauss_vels(i,:))]};
      endif
    end
  end

  %Plot ZAMS, only one of them
  sub_folder = strcat(gauss_fields(1,:), '_', rotational_vels(1,:));
  full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
  zams = calculate_ZAMS(full_path);
  line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 3, "linestyle", "--", "color", "k");

  % Plot sun reference
  plot(sun_age, log10(sun_T_eff), '*', 'markersize', 15, 'color', [0.5,0.1,0.8]);
  plot(sun_age, sun_log_g, 's', 'markersize', 10, 'color', [0.5,0.1,0.8], 'markerfacecolor', [0.5,0.1,0.8]);

  % Plot reference marks
  line("xdata",[2.5e7,2.5e7], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 2, "linestyle", ":", "color", "cyan");
  line("xdata",[3.5e7,3.5e7], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 2, "linestyle", ":", "color", "cyan");
  line("xdata",[5.4e7,5.4e7], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 2, "linestyle", ":", "color", "magenta");
  line("xdata",[11.2e7,11.2e7], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 2, "linestyle", ":", "color", "magenta");

  l = legend(labels, "location", leg_loc);
  set (l, "fontsize", legend_font_size);
  %legend boxoff
  xlabel('star age (Ga)', 'fontsize', axis_font_size);
  ylabel('T_{eff} & log g', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');

  plot_path = strcat(tables_parent_folder, '/', aidx);
  save_figure(f, strcat(plot_path, '/', afilename, aidx));
end


function age_vs_radius_plots(gauss_fields, rotational_vels, is_var_vel, ytick, axis_limits, leg_loc, atitle, afilename, aidx)
  global data_parent_folder;
  global tables_parent_folder;
  global filename;
  global star_age_col;
  global header_lines;
  global log_R_col
  global sun_age;
  global sun_radius;
  global sun_vel_rot;
  global colors;
  global sun_T_eff;
  global sun_log_g;
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

      fmt = get_parsing_fmt([star_age_col, log_R_col]);

      A = read_matrix_from_file(full_path, fmt, header_lines, 2);

      plot_radius(A, colors(i*j,:), line_width, ytick, axis_limits);


      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['log R -', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['log R -', strtrim(gauss_fields(i,:))]};
      endif
    end
  end

  %Plot ZAMS, only one of them
  sub_folder = strcat(gauss_fields(1,:), '_', rotational_vels(1,:));
  full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
  zams = calculate_ZAMS(full_path);
  line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 3, "linestyle", "--", "color", "k");

  % Plot sun reference
  plot(sun_age, 0, '*', 'markersize', 15, 'color', [0.5,0.1,0.8]);

  l = legend(labels, "location", leg_loc);
  set (l, "fontsize", legend_font_size);
  %legend boxoff
  xlabel('star age (Ga)', 'fontsize', axis_font_size);
  ylabel('log (R/R_{sun})', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');

  plot_path = strcat(tables_parent_folder, '/', aidx);
  save_figure(f, strcat(plot_path, '/', afilename, aidx));
end




function age_vs_cz_size_plots_special(rotational_vels, is_var_vel, ytick, axis_limits, leg_loc, atitle, afilename)
  global data_parent_folder;
  global filename;
  global star_age_col;
  global log_R_col;
  global sz_top_radius_col;
  global sz_bot_radius_col;
  global header_lines;
  global sun_age;
  global sun_vel_rot;
  %global colors;
  global title_font_size;
  global axis_font_size;
  global legend_font_size;
  global line_width;

  colors = ['k'; 'r'; 'g'; 'b'; 'y'; 'm'; 'k'; 'r'; 'g'];

  hold('on');
  labels = {};

  f = format_figure();

  gauss_fields = ['0g'; '3g'; '3.5g'; '4g'; '4.5g'; '5g'; '5gx2s'; '5gx2st'; '5gx4s'];

  for i=1:rows(gauss_fields)
    for j=1:rows(rotational_vels)
      sub_folder = strcat(gauss_fields(i,:), '_', rotational_vels(j,:));
      full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);

      fmt = get_parsing_fmt([star_age_col, log_R_col, sz_top_radius_col, sz_bot_radius_col]);

      A = read_matrix_from_file(full_path, fmt, header_lines, 4);

      if (i == 7)
        plot_size_cz2(A, colors(i*j,:), line_width, ytick, axis_limits, -0.005, ':');
      elseif (i == 8)
        plot_size_cz2(A, colors(i*j,:), line_width, ytick, axis_limits, -0.0065, ':');
      elseif (i == 9)
        plot_size_cz2(A, colors(i*j,:), line_width, ytick, axis_limits, -0.0090, ':');
      else
        plot_size_cz2(A, colors(i*j,:), line_width, ytick, axis_limits, 0, '-');
      endif

      % Plot ZAMS reference
      %zams = calculate_ZAMS(full_path);
      %line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 2, "linestyle", "--", "color", colors(i*j,:))

      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['CZ_{rad}-', strtrim(rotational_vels(j,:))]};
        %labels = {labels{:}, ['ZAMS-', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['CZ_{rad}-', strtrim(gauss_fields(i,:))]};
        %labels = {labels{:}, ['ZAMS-', strtrim(gauss_fields(i,:))]};
      endif
    end
  end

  %Plot ZAMS, only one of them
  sub_folder = strcat(gauss_fields(1,:), '_', rotational_vels(1,:));
  full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
  zams = calculate_ZAMS(full_path);
  line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 3, "linestyle", "--", "color", "k");

  l = legend(labels, "location", leg_loc);
  set (l, "fontsize", legend_font_size);
  %legend boxoff
  xlabel('star age (Ga)', 'fontsize', axis_font_size);
  ylabel('Size conv. zone (radius conv. zone/radius star)', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');

  save_figure(f, afilename);
end


function age_vs_m_dot_plots(gauss_fields, rotational_vels, is_var_vel, ytick, axis_limits, leg_loc, atitle, afilename, aidx)
  global data_parent_folder;
  global tables_parent_folder;
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
      %zams = calculate_ZAMS(full_path);
      %line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 2, "linestyle", "--", "color", colors(i*j,:))

      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['Mass loss-', strtrim(rotational_vels(j,:))]};
        %labels = {labels{:}, ['ZAMS-', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['Mass loss-', strtrim(gauss_fields(i,:))]};
        %labels = {labels{:}, ['ZAMS-', strtrim(gauss_fields(i,:))]};
      endif

    end
  end

  %Plot ZAMS, only one of them
  sub_folder = strcat(gauss_fields(1,:), '_', rotational_vels(1,:));
  full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
  zams = calculate_ZAMS(full_path);
  line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 3, "linestyle", "--", "color", "k");

  l = legend(labels, "location", leg_loc);
  set (l, "fontsize", legend_font_size);
  %legend boxoff
  xlabel('star age (Ga)', 'fontsize', axis_font_size);
  ylabel('Mass loss (log(M_{sun}/yr))', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');

  plot_path = strcat(tables_parent_folder, '/', aidx);
  save_figure(f, strcat(plot_path, '/', afilename, aidx));
end

function age_vs_reimers_m_dot_plots(gauss_fields, rotational_vels, is_var_vel, ytick, axis_limits, atitle, afilename)
  global data_parent_folder;
  global filename;
  global star_age_col;
  global star_mass_col;
  global log_L_col;
  global log_R_col;
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

      fmt = get_parsing_fmt([star_age_col, star_mass_col, log_L_col, log_R_col]);

      A = read_matrix_from_file(full_path, fmt, header_lines, 4);

      plot_reimers_m_dot(A, colors(i*j,:), line_width, ytick, axis_limits);

      % Plot ZAMS reference
      %zams = calculate_ZAMS(full_path);
      %line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 2, "linestyle", "--", "color", colors(i*j,:))

      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['Mass loss-', strtrim(rotational_vels(j,:))]};
        %labels = {labels{:}, ['ZAMS-', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['Mass loss-', strtrim(gauss_fields(i,:))]};
        %labels = {labels{:}, ['ZAMS-', strtrim(gauss_fields(i,:))]};
      endif

    end
  end

  %Plot ZAMS, only one of them
  sub_folder = strcat(gauss_fields(1,:), '_', rotational_vels(1,:));
  full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
  zams = calculate_ZAMS(full_path);
  line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 3, "linestyle", "--", "color", "k");

  l = legend(labels, "location", "southeastoutside");
  set (l, "fontsize", legend_font_size);
  legend boxoff
  xlabel('star age (yrs)', 'fontsize', axis_font_size);
  ylabel('Mass loss (log(solar mass/yr))', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');

  save_figure(f, afilename);
end


function age_vs_m_l_r_plots(gauss_fields, rotational_vels, is_var_vel, ytick, x_limits, atitle, afilename)
    global data_parent_folder;
  global filename;
  global star_age_col;
  global star_mass_col;
  global log_L_col;
  global log_R_col;
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

      fmt = get_parsing_fmt([star_age_col, star_mass_col, log_L_col, log_R_col]);

      A = read_matrix_from_file(full_path, fmt, header_lines, 4);
      %Get the index of the last records lower than or equal to the temporal limits
      ix_ini = find(A(:,1)<=x_limits(1), 1, 'last');
      ix_end = find(A(:,1)<=x_limits(2), 1, 'last');

      %Calculate maximum for y axis
      %Get vel max value, divide it by 10, plus 1, multiply by 10
      ymax = (idivide(power(10, max(A(ix_ini:ix_end,3))), ytick, "fix") + 1) * ytick;

      %Get vel max value, divide it by 10, minus 1, multiply by 10
      ymin = (idivide(power(10, min(A(ix_ini:ix_end,3))), ytick, "fix") - 1) * ytick;

      plot_m_l_r(A, colors(i*j,:), line_width, ytick, [int64(x_limits(1)), int64(x_limits(2)), ymin, ymax]);
      %plot_m_l_r(A, colors(i*j,:), line_width, ytick, x_limits);

      % Plot ZAMS reference
      %zams = calculate_ZAMS(full_path);
      % Here we can properly assing the ymin and ymax values for all the plots
      % Each of them will potentially have a different one. We opt for fixing
      % the limits by hand.
      %line("xdata",[zams,zams], "ydata",[-10,200], "linewidth", 2, "linestyle", "--", "color", colors(i*j,:))


      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['Mass-', strtrim(rotational_vels(j,:))]};
        labels = {labels{:}, ['Luminosity-', strtrim(rotational_vels(j,:))]};
        labels = {labels{:}, ['Radius-', strtrim(rotational_vels(j,:))]};
        %labels = {labels{:}, ['ZAMS-', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['Mass-', strtrim(gauss_fields(i,:))]};
        labels = {labels{:}, ['Luminosity-', strtrim(gauss_fields(i,:))]};
        labels = {labels{:}, ['Radius-', strtrim(gauss_fields(i,:))]};
        %labels = {labels{:}, ['ZAMS-', strtrim(gauss_fields(i,:))]};
      endif

    end
  end

  %Plot ZAMS, only one of them
  sub_folder = strcat(gauss_fields(1,:), '_', rotational_vels(1,:));
  full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
  zams = calculate_ZAMS(full_path);
  line("xdata",[zams,zams], "ydata",[-10,200], "linewidth", 3, "linestyle", "--", "color", "k");

  l = legend(labels, "location", "southeastoutside");
  %l = legend(labels, "location", "southoutside", "orientation", "horizontal");
  set (l, "fontsize", legend_font_size);
  legend boxoff
  xlabel('star age (yrs)', 'fontsize', axis_font_size);
  ylabel('M, L, R (solar units)', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');
  save_figure(f, afilename);

end



function age_vs_vel_plots(gauss_fields, rotational_vels, is_var_vel, ytick, x_limits, leg_loc, atitle, afilename, aidx)
  global tables_parent_folder;
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
      ymax = (idivide(max(A(ix_ini:ix_end,2)), int16(ytick), "fix") + 1) * ytick;

      %Get vel max value, divide it by 10, minus 1, multiply by 10
      ymin = (idivide(min(A(ix_ini:ix_end,2)), int16(ytick), "fix") - 1) * ytick;

      plot_vel_rot(A, colors(i*j,:), line_width, ytick, [int64(x_limits(1)), int64(x_limits(2)), ymin, ymax]);

      % Plot ZAMS reference
      %zams = calculate_ZAMS(full_path);
      % Here we can properly assing the ymin and ymax values for all the plots
      % Each of them will potentially have a different one. We opt for fixing
      % the limits by hand.
      %line("xdata",[zams,zams], "ydata",[-10,200], "linewidth", 2, "linestyle", "--", "color", colors(i*j,:));


      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['Surface-', strtrim(rotational_vels(j,:))]};
        %Dont show lim sup CZ
        %labels = {labels{:}, ['Top CZ-', strtrim(rotational_vels(j,:))]};
        labels = {labels{:}, ['Bottom CZ-', strtrim(rotational_vels(j,:))]};
        %Dont show individual ZAMS
        %labels = {labels{:}, ['ZAMS-', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['Surface-', strtrim(gauss_fields(i,:))]};
        %Dont show lim sup CZ
        %labels = {labels{:}, ['Top CZ-', strtrim(gauss_fields(i,:))]};
        labels = {labels{:}, ['Bottom CZ-', strtrim(gauss_fields(i,:))]};
        %Dont show individual ZAMS
        %labels = {labels{:}, ['ZAMS-', strtrim(gauss_fields(i,:))]};
      endif

    end
  end
  %Plot ZAMS, only one of them
  sub_folder = strcat(gauss_fields(1,:), '_', rotational_vels(1,:));
  full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
  zams = calculate_ZAMS(full_path);
  % Here we can properly assing the ymin and ymax values for all the plots
  % Each of them will potentially have a different one. We opt for fixing
  % the limits by hand.
  line("xdata",[zams,zams], "ydata",[-10,500], "linewidth", 3, "linestyle", "--", "color", "k");


  % Plot sun reference
  plot(sun_age, sun_vel_rot, '*', 'markersize', 15, 'linewidth', 2, 'color', [0.5,0.1,0.8]);

  l = legend(labels, "location", leg_loc);

  set (l, "fontsize", legend_font_size);
  %legend boxoff
  xlabel('star age (Ga)', 'fontsize', axis_font_size);
  ylabel('Rotational Vel (km/s)', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');

  plot_path = strcat(tables_parent_folder, '/', aidx);
  save_figure(f, strcat(plot_path, '/', afilename, aidx));

end


function age_vs_omega_plots(gauss_fields, rotational_vels, is_var_vel, ytick, x_limits, leg_loc, atitle, afilename, aidx)
  global tables_parent_folder;
  global data_parent_folder;
  global filename;
  global star_age_col;
  global surf_avg_omega_col;
  global sz_top_omega_col;
  global sz_bot_omega_col;
  global core_top_omega_col;
  global header_lines;
  global sun_age;
  global sun_omega;
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

      fmt = get_parsing_fmt([star_age_col, surf_avg_omega_col, sz_top_omega_col, sz_bot_omega_col, core_top_omega_col]);

      A = read_matrix_from_file(full_path, fmt, header_lines, 5);
      %Get the index of the last records lower than or equal to the temporal limits
      ix_ini = find(A(:,1)<=x_limits(1), 1, 'last');
      ix_end = find(A(:,1)<=x_limits(2), 1, 'last');

      %Calculate maximum for y axis
      %Get vel max value, divide it by ytick, plus 1, multiply by ytick
      ymax = (idivide(max(A(ix_ini:ix_end,2)), int8(ytick), "fix") + 1) * ytick;
      %ymax = ymax / sun_omega;
      ymax = 400.0;

      %Get vel max value, divide it by ytick, minus 1, multiply by ytick
      ymin = (idivide(min(A(ix_ini:ix_end,2)), int8(ytick), "fix") - 1) * ytick;
      ymin = ymin / sun_omega;
      %ymin = 0.5;

      plot_omega(A, colors(i*j,:), line_width, ytick, [int64(x_limits(1)), int64(x_limits(2)), ymin, ymax]);

      % Plot ZAMS reference
      %zams = calculate_ZAMS(full_path);
      % Here we can properly assing the ymin and ymax values for all the plots
      % Each of them will potentially have a different one. We opt for fixing
      % the limits by hand.
      %line("xdata",[zams,zams], "ydata",[-10,200], "linewidth", 2, "linestyle", "--", "color", colors(i*j,:));


      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['Surface-', strtrim(rotational_vels(j,:))]};
        %Dont show lim sup CZ
        %labels = {labels{:}, ['Top CZ-', strtrim(rotational_vels(j,:))]};
        labels = {labels{:}, ['Top core-', strtrim(rotational_vels(j,:))]};
        %Dont show individual ZAMS
        %labels = {labels{:}, ['ZAMS-', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['Surface-', strtrim(gauss_fields(i,:))]};
        %Dont show lim sup CZ
        %labels = {labels{:}, ['Top CZ-', strtrim(gauss_fields(i,:))]};
        labels = {labels{:}, ['Top core-', strtrim(gauss_fields(i,:))]};
        %Dont show individual ZAMS
        %labels = {labels{:}, ['ZAMS-', strtrim(gauss_fields(i,:))]};
      endif

    end
  end
  %Plot ZAMS, only one of them
  sub_folder = strcat(gauss_fields(1,:), '_', rotational_vels(1,:));
  full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
  zams = calculate_ZAMS(full_path);

  % Here we can properly assing the ymin and ymax values for all the plots
  % Each of them will potentially have a different one. We opt for fixing
  % the limits by hand.
  line("xdata",[zams,zams], "ydata",[0.001,100], "linewidth", 3, "linestyle", "--", "color", "k");


  % Plot sun reference
  plot(sun_age, sun_omega, '*', 'markersize', 15, 'linewidth', 2, 'color', [0.5,0.1,0.8]);

  l = legend(labels, "location", leg_loc);

  set (l, "fontsize", legend_font_size);
  %legend boxoff
  xlabel('star age (Ga)', 'fontsize', axis_font_size);
  ylabel('Omega star/Omega Sun', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');
  plot_path = strcat(tables_parent_folder, '/', aidx);
  save_figure(f, strcat(plot_path, '/', afilename, aidx));
end

function age_vs_alpha_mlt(gauss_fields, rotational_vels, is_var_vel, ytick, axis_limits, leg_loc, atitle, afilename, aidx)
  global data_parent_folder;
  global tables_parent_folder;
  global filename;
  global star_age_col;
  global alpha_mlt_col;
  global header_lines;
  global colors;
  global title_font_size;
  global axis_font_size;
  global legend_font_size;
  global line_width;
  global sun_age;
  amlt_samadi_2006 = 1.76;
  amlt_sonoi_2019 = 1.78;

  hold('on');
  labels = {};

  f = format_figure();

  for i=1:rows(gauss_fields)
    for j=1:rows(rotational_vels)
      sub_folder = strcat(gauss_fields(i,:), '_', rotational_vels(j,:));
      full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);

      fmt = get_parsing_fmt([star_age_col, alpha_mlt_col]);

      A = read_matrix_from_file(full_path, fmt, header_lines, 2);
      plot_alpha_mlt(A, colors(i*j,:), line_width, ytick, axis_limits);


      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['\alpha_{MLT}-', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['\alpha_{MLT}-', strtrim(gauss_fields(i,:))]};
      endif

    end
  end
  %Plot ZAMS, only one of them
  sub_folder = strcat(gauss_fields(1,:), '_', rotational_vels(1,:));
  full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
  zams = calculate_ZAMS(full_path);
  line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 3, "linestyle", "--", "color", "k");

  %Plot amlt
  plot(sun_age, amlt_samadi_2006, 'd', 'markersize', 15, 'color', [0.5,0.1,0.8]);
  plot(sun_age, amlt_sonoi_2019, '*', 'markersize', 15, 'color', [0.5,0.1,0.8]);

  l = legend(labels, "location", leg_loc);

  set (l, "fontsize", legend_font_size);
  set(l, 'interpreter', 'tex');

  %legend boxoff
  xlabel('star age (Ga)', 'fontsize', axis_font_size);
  ylabel('\alpha_{MLT}', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');

  plot_path = strcat(tables_parent_folder, '/', aidx);
  plot_path
  save_figure(f, strcat(plot_path, '/', afilename, aidx));

end

function omegs_vs_mag_field(gauss_fields, rotational_vels, show_limits, ytick, axis_limits, leg_loc, atitle, afilename, aidx)
  global data_parent_folder;
  global tables_parent_folder;
  global filename;
  global star_age_col;
  global surf_avg_omega_col;
  global tau_c_col;
  global rossby_col;
  global rossby_norm_col;
  global b_equi_col;
  global b_star_col;
  global f_min_col;
  global f_max_col;
  global f_star_col;
  global bf_min_col;
  global bf_max_col;
  global bf_star_col;
  global alfven_r_col;
  global sun_vel_rot;
  global sun_gauss_field;
  global sun_age;

  global header_lines;
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

      fmt = get_parsing_fmt([star_age_col, surf_avg_omega_col, bf_min_col, bf_max_col, bf_star_col]);

      A = read_matrix_from_file(full_path, fmt, header_lines, 5);

      plot_mag_field(A, colors(i*j,:), line_width, show_limits, ytick, axis_limits);


      %Generate serie labels

      labels = {labels{:}, ['\omega_{0} -', strtrim(rotational_vels(j,:))]};
      if (show_limits)
        labels = {labels{:}, ['B_{f_{min}} -', strtrim(rotational_vels(j,:))]};
        labels = {labels{:}, ['B_{f_{max}} -', strtrim(rotational_vels(j,:))]};
      endif
      labels = {labels{:}, ['B_{f_{mean}} -', strtrim(rotational_vels(j,:))]};
    end
  end
  %Plot ZAMS, only one of them
  sub_folder = strcat(gauss_fields(1,:), '_', rotational_vels(1,:));
  full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
  zams = calculate_ZAMS(full_path);
  line("xdata",[zams,zams], "ydata",[axis_limits(3),axis_limits(4)], "linewidth", 3, "linestyle", "--", "color", "k");

  % Plot sun reference
  %plot(sun_age, sun_vel_rot, '*', 'markersize', 15, 'linewidth', 2, 'color', [0.5,0.1,0.8]);
  %plot(sun_age, sun_gauss_field, 's', 'markersize', 10, 'color', [0.5,0.1,0.8], 'markerfacecolor', [0.5,0.1,0.8]);

  l = legend(labels, "location", leg_loc);

  set(l, "fontsize", legend_font_size);
  set(l, 'interpreter', 'tex');
  %legend boxoff
  xlabel('star age (Ga)', 'fontsize', axis_font_size);
  ylabel('Mean mag. field (G) & \Omega/\Omega_{Sun}', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');

  plot_path = strcat(tables_parent_folder, '/', aidx);
  save_figure(f, strcat(plot_path, '/', afilename, aidx));
end





function kipperhahn_plots(gauss_fields, rotational_vels, is_var_vel, ytick, x_limits, leg_loc, atitle, afilename)
  global data_parent_folder;
  global filename;
  global star_age_col;
  global header_lines;
  global colors;
  global title_font_size;
  global axis_font_size;
  global legend_font_size;
  global line_width;
  global num_mix_regions;
  global mix_type_ini_col;
  global mix_relr_top_ini_col;


  hold('on');
  labels = {};

  f = format_figure();

      sub_folder = strcat(gauss_fields(1,:), '_', rotational_vels(1,:));
      full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);

      mix_type_cols = [star_age_col];
      mix_size_cols = [star_age_col];
      for k=0:num_mix_regions-1
        mix_type_cols(end+1) = mix_type_ini_col + (k*2);
        mix_size_cols(end+1) = mix_relr_top_ini_col + (k*2);
      end
      mix_type_cols
      fmt = get_parsing_fmt(mix_type_cols);
      A = read_matrix_from_file(full_path, fmt, header_lines, num_mix_regions+1);

      fmt2 = get_parsing_fmt(mix_size_cols);
      B = read_matrix_from_file(full_path, fmt2, header_lines, num_mix_regions+1);

      ymax = 1.0;
      ymin = 0.0;

      plot_kipperhahn(A, B, 1, line_width, ytick, [int64(x_limits(1)), int64(x_limits(2)), ymin, ymax]);

      %plot_omega(A, colors(i*j,:), line_width, ytick, [x_limits(1), x_limits(2), ymin, ymax]);
      %line("xdata",[zams,zams], "ydata",[-10,200], "linewidth", 2, "linestyle", "--", "color", colors(i*j,:));


      %Generate serie labels
      if (is_var_vel)
        %labels = {labels{:}, ['Surface-', strtrim(rotational_vels(1,:))]};
        %Dont show lim sup CZ
        %labels = {labels{:}, ['Top CZ-', strtrim(rotational_vels(j,:))]};
        %labels = {labels{:}, ['Top core-', strtrim(rotational_vels(1,:))]};
        %Dont show individual ZAMS
        %labels = {labels{:}, ['ZAMS-', strtrim(rotational_vels(j,:))]};
      else
        %labels = {labels{:}, ['Surface-', strtrim(gauss_fields(1,:))]};
        %Dont show lim sup CZ
        %labels = {labels{:}, ['Top CZ-', strtrim(gauss_fields(i,:))]};
        %labels = {labels{:}, ['Top core-', strtrim(gauss_fields(1,:))]};
        %Dont show individual ZAMS
        %labels = {labels{:}, ['ZAMS-', strtrim(gauss_fields(i,:))]};
      endif

  %Plot ZAMS, only one of them
  sub_folder = strcat(gauss_fields(1,:), '_', rotational_vels(1,:));
  full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
  zams = calculate_ZAMS(full_path);
  % Here we can properly assing the ymin and ymax values for all the plots
  % Each of them will potentially have a different one. We opt for fixing
  % the limits by hand.
  line("xdata",[zams,zams], "ydata",[0.0,1.0], "linewidth", 5, "linestyle", "--", "color", "g");

  %l = legend(labels, "location", leg_loc);
  legend off;

  %set (l, "fontsize", legend_font_size);
  %legend boxoff
  xlabel('star age (yrs)', 'fontsize', axis_font_size);
  ylabel('Mixing type', 'fontsize', axis_font_size);
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
  %print(f,'-dpdflatexstandalone','-color',[title,'.pdf']);
  print(f,'-deps','-color',[title,'.eps']);
  close;
  %print(f,'-dpng','-color',[title,'.png']);
  %$ for i in *.eps; do ps2pdf $i $(basename -s .eps $i).pdf ; done
end

function plot_mb_activation(A, color, width, axis_limits)
  global tick_font_size;

  %Plot values
  plot(A(:,1), A(:,3), color, 'linewidth', width);
  plot(A(:,1), A(:,2), color, 'linewidth', width, 'linestyle', '--');

  %Axis scales
  axis_limits
  %set(gca, 'XScale', 'log', 'XTick', 1.0e0:1.0e1:1.0e10);
  set(gca, 'XScale', 'log');
  grid on;
  set(gca,'xminorgrid','off');

  %Axis limits
  set(gca,'YTick',0.0:1:axis_limits(4));

  %Axis ticks
  axis(axis_limits);
  xticks = get (gca, "xtick");
  %xlabels = arrayfun (@(x) sprintf ("%.e", x), xticks, "uniformoutput", false);
  %xlabels = [' ';'1.00e+02';' ';'1.00e+04';' ';'1.00e+06';' ';'1.00e+08';' ';'10.00e+10'];
  xlabels = [' ';'0.0000001';' ';'0.00001';' ';'0.001';' ';'0.1';' ';'10.0'];
  set (gca, "xticklabel", xlabels) ;


  yticks = get (gca, "ytick");
  ylabels = arrayfun (@(x) sprintf ("%2.1f", x), yticks, "uniformoutput", false);
  set (gca, "yticklabel", ylabels);
  set(gca, 'fontsize', tick_font_size);

end


function age_vs_mb_activation(gauss_fields, rotational_vels, is_var_vel, leg_loc, atitle, afilename, aidx)
  global data_parent_folder;
  global tables_parent_folder;
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
      %zams = calculate_ZAMS(full_path);
      % Here we can properly assing the ymin and ymax values for all the plots
      % Each of them will potentially have a different one. We opt for fixing
      % the limits by hand.
      %line("xdata",[zams,zams], "ydata",[-10,200], "linewidth", 2, "linestyle", "--", "color", colors(i*j,:))


      %Generate serie labels
      if (is_var_vel)
        labels = {labels{:}, ['MB-', strtrim(rotational_vels(j,:))]};
        labels = {labels{:}, ['CR-', strtrim(rotational_vels(j,:))]};
        %labels = {labels{:}, ['ZAMS-', strtrim(rotational_vels(j,:))]};
      else
        labels = {labels{:}, ['MB-', strtrim(gauss_fields(i,:))]};
        labels = {labels{:}, ['CR-', strtrim(gauss_fields(i,:))]};
        %labels = {labels{:}, ['ZAMS-', strtrim(gauss_fields(i,:))]};
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

  %Plot ZAMS, only one of them
  sub_folder = strcat(gauss_fields(1,:), '_', rotational_vels(1,:));
  full_path = strcat(data_parent_folder, '/', sub_folder, '/', filename);
  zams = calculate_ZAMS(full_path);
  line("xdata",[zams,zams], "ydata",[-10,200], "linewidth", 3, "linestyle", "--", "color", "k");



  yticks(y_ticks);
  yticklabels(y_labels);

  l = legend(labels, "location", leg_loc);
  set (l, "fontsize", legend_font_size);
  legend boxoff
  xlabel('star age (Ga)', 'fontsize', axis_font_size);
  ylabel('MB activation & Radiative vs. Convective core', 'fontsize', axis_font_size);
  title(atitle, 'fontsize', title_font_size);

  hold('off');

  plot_path = strcat(tables_parent_folder, '/', aidx);
  save_figure(f, strcat(plot_path, '/', afilename, aidx));
end



function plot_0G_var_vel(rot_vels,idx)
  global gauss_fields;
  global idx_0_0G;

  age_vs_li_plots(gauss_fields(idx_0_0G,:), rot_vels, true, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', '', 'li_var_vel_0_0g_', num2str(idx));
end


function plot_0G_var_vel_z1(rot_vels,idx)
  global gauss_fields;
  global idx_0_0G;

  age_vs_li_plots(gauss_fields(idx_0_0G,:), rot_vels, true, 0.1, [1.0e7, 1.0e8, 2.1, 2.4], 'northeast', 'A(Li7) - 0G & var. rotational velocity', 'li_var_vel_0_0g_z1_', num2str(idx));
end

function plot_2G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_2_0G;
  gauss_fields = ['2g'];
  %rotational_vels = ['037crit';'0375crit';'038crit';'0385crit';'039crit';'0395crit'];
  %rotational_vels = ['035crit';'034crit';'033crit';'032crit';'031crit'];
  rotational_vels = ['026crit';'027crit';'028crit';'029crit';'030crit'];

  age_vs_li_plots(gauss_fields(idx_2_0G,:), rot_vels, true, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', 'A(Li7) - 2.0G & var. rotational velocity', 'li_var_vel_2_0g', num2str(idx));
end

function plot_2_5G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_2_5G;

  age_vs_li_plots(gauss_fields(idx_2_5G,:), rot_vels, true, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', 'A(Li7) - 2.5G & var. rotational velocity', 'li_var_vel_2_5g', num2str(idx));
end


function plot_3_0G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_3_0G;

  age_vs_li_plots(gauss_fields(idx_3_0G,:), rot_vels, true, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', 'A(Li7) - 3.0G & var. rotational velocity', 'li_var_vel_3_0g', num2str(idx));
end

function plot_3_0G_0314vc(rot_vels, idx)
  global gauss_fields;
  global idx_3_0G;

  age_vs_li_plots(gauss_fields(idx_3_0G,:), rot_vels, true, 0.5, [1.0e5,1.0e10,0,3.5], 'west', 'A(Li7) - 3.0G & vcrit=0.0314', 'li_3_0g_0314vc', num2str(idx));
end



function plot_3_3G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_3_3G;

  age_vs_li_plots(gauss_fields(idx_3_3G,:), rot_vels, true, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', 'A(Li7) - 3.3G & var. rotational velocity', 'li_var_vel_3_3g', num2str(idx));
end


function plot_3_5G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_3_5G;
  %rotational_vels2 = ['029crit';'031crit';'033crit'];

  age_vs_li_plots(gauss_fields(idx_3_5G,:), rot_vels, true, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', 'A(Li7) - 3.5G & var. rotational velocity', 'li_var_vel_3_5g', num2str(idx));
end

function plot_4_0G_var_vel(rot_vels,idx)
  global gauss_fields;
  global idx_4_0G;

  age_vs_li_plots(gauss_fields(idx_4_0G,:), rot_vels, true, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', '', 'li_var_vel_4_0g', num2str(idx));
end

function plot_4_0G_var_vel_st(rot_vels,idx)
  global gauss_fields;
  global idx_4_0G;
  rotational_vels = ['0196crit'; '0196crit_st'];
  age_vs_li_plots(gauss_fields(idx_4_0G,:), rot_vels, true, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', '', 'li_var_vel_4_0g', num2str(idx));
end

function plot_4_0G_var_vel_vc5(rot_vels, idx)
  global gauss_fields;
  rotational_vels = ['0196crit'; '0196crit_vc5'; '0196crit_vc5_md5'];
  age_vs_li_plots(gauss_fields(idx_4_0G,:), rot_vels, true, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', 'A(Li7) - 4.0G & var. rotational velocity', 'li_var_vel_4_0g', num2str(idx));
end

function plot_4_3G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_4_3G;

  age_vs_li_plots(gauss_fields(idx_4_3G,:), rot_vels, true, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', 'A(Li7) - 4.3G & var. rotational velocity', 'li_var_vel_4_3g', num2str(idx));
end


function plot_4_5G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_4_5G;

  age_vs_li_plots(gauss_fields(idx_4_5G,:), rot_vels, true, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', 'A(Li7) - 4.5G & var. rotational velocity', 'li_var_vel_4_5g', num2str(idx));
end

function plot_4_5G_var_vel_special(rot_vels, idx)
  global gauss_fields;
  global idx_4_5G;

  age_vs_li_plots(gauss_fields(idx_4_5G,:), rot_vels, true, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', 'A(Li7) - 4.5G & var. rotational velocity', 'li_var_vel_4_5g', num2str(idx));
end


function plot_5_0G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_5_0G;

  age_vs_li_plots(gauss_fields(idx_5_0G,:), rot_vels, true, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', 'A(Li7) - 5.0G & var. rotational velocity', 'li_var_vel_5_0g', num2str(idx));
end

function plot_5_5G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_5_5G;

  age_vs_li_plots(gauss_fields(idx_5_5G,:), rot_vels, true, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', 'A(Li7) - 5.5G & var. rotational velocity', 'li_var_vel_5_5g', num2str(idx));
end

function plot_XG_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_X_G;

  age_vs_li_plots(gauss_fields(idx_X_G,:), rot_vels, true, 0.5, [1.0e5,1.0e10,0.8,3.5], 'southwest', '', 'li_var_vel_var_g_', num2str(idx));
end



function plot_0084vc_var_g(mag_fields, idx)
  global rotational_vels;
  global idx_0084crit;

  age_vs_li_plots(mag_fields, rotational_vels(idx_0084crit,:), false, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', 'A(Li7) - vcrit=0.0084 & var. magnetic field', 'li_vc_0084_var_g', num2str(idx));
end

function plot_014vc_var_g(mag_fields, idx)
  global rotational_vels;
  global idx_014crit;

  age_vs_li_plots(mag_fields, rotational_vels(idx_014crit,:), false, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', 'A(Li7) - vcrit=0.014 & var. magnetic field', 'li_vc_014_var_g', num2str(idx));
end
function plot_0196vc_var_g(mag_fields, idx)
  global rotational_vels;
  global idx_0196crit;

  age_vs_li_plots(mag_fields, rotational_vels(idx_0196crit,:), false, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', 'A(Li7) - vcrit=0.0196 & var. magnetic field', 'li_vc_0196_var_g', num2str(idx));
end
function plot_028vc_var_g(mag_fields, idx)
  global rotational_vels;
  global idx_028crit;

  age_vs_li_plots(mag_fields, rotational_vels(idx_028crit,:), false, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', 'A(Li7) - vcrit=0.028 & var. magnetic field', 'li_vc_028_var_g', num2str(idx));
end

function plot_0336vc_var_g(mag_fields, idx)
  global rotational_vels;
  global idx_0336crit;

  age_vs_li_plots(mag_fields, rotational_vels(idx_0336crit,:), false, 0.5, [1.0e5,1.0e10,0,3.5], 'southwest', 'A(Li7) - vcrit=0.0336 & var. magnetic field', 'li_vc_0336_var_g', num2str(idx));
end



function plot_vel_rot_0G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_0_0G;

  age_vs_vel_plots(gauss_fields(idx_0_0G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northwest', '', 'rot_vel_var_vel_0_0g_', num2str(idx));
end

function plot_vel_rot_2G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_2_0G;

  age_vs_vel_plots(gauss_fields(idx_2_0G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northwest', 'Rotational vel - 2.0G & var. rotational velocity', 'rot_vel_var_vel_2_0g', num2str(idx));
end

function plot_vel_rot_2_5G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_2_5G;

  age_vs_vel_plots(gauss_fields(idx_2_5G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northwest', 'Rotational vel - 2.5G & var. rotational velocity', 'rot_vel_var_vel_2_5g', num2str(idx));
end


function plot_vel_rot_3_0G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_3_0G;

  age_vs_vel_plots(gauss_fields(idx_3_0G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northwest', 'Rotational vel - 3.0G & var. rotational velocity', 'rot_vel_var_vel_3_0g', num2str(idx));
end

function plot_vel_rot_3_0G_0314vc(rot_vels, idx)
  global gauss_fields;
  global idx_3_0G;

  age_vs_vel_plots(gauss_fields(idx_3_0G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northwest', 'Rotational vel - 3.0G & vcrit=0.0314', 'rot_vel_3_0g_0314vc', num2str(idx));
end

function plot_vel_rot_3_0G_var_vel_mlt(rot_vels, idx)
  global gauss_fields;
  global idx_3_0G;

  age_vs_vel_plots(gauss_fields(idx_3_0G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northwest', 'Rotational vel - 3.0G & var. rotational velocity', 'rot_vel_var_vel_3_0g', num2str(idx));
end


function plot_vel_rot_3_3G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_3_3G;

  age_vs_vel_plots(gauss_fields(idx_3_3G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northwest', 'Rotational vel - 3.3G & var. rotational velocity', 'rot_vel_var_vel_0_0g', num2str(idx));
end



function plot_vel_rot_3_5G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_3_5G;

  age_vs_vel_plots(gauss_fields(idx_3_5G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northwest', 'Rotational velocity - 3.5G & var. rotational velocity', 'rot_vel_var_vel_3_5g', num2str(idx));
end

function plot_vel_rot_4G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_4_0G;

  age_vs_vel_plots(gauss_fields(idx_4_0G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northwest', '', 'rot_vel_var_vel_4_0g', num2str(idx));
end

function plot_vel_rot_4G_var_vel_vc5_md5(rot_vels, idx)
  global gauss_fields;
  global idx_4_0G;
  rotational_vels = ['0196crit'; '0196crit_vc5'; '0196crit_vc5_md5'];

  age_vs_vel_plots(gauss_fields(idx_4_0G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northwest', 'Rotational velocity - 4.0G & var. rotational velocity', 'rot_vel_var_vel_4_0g', num2str(idx));
end


function plot_vel_rot_4G_var_vel_z1(rot_vels, idx)
  global gauss_fields;
  global idx_4_0G;

  age_vs_vel_plots(gauss_fields(idx_4_0G,:), rot_vels, true, 2, [3.0e9, 1.0e10], 'northeast', '', 'rot_vel_var_vel_4_0g_z1', num2str(idx));
end

function plot_vel_rot_4G_var_vel_z1_vc5_md5(rot_vels, idx)
  global gauss_fields;
  global idx_4_0G;
  rotational_vels = ['0196crit'; '0196crit_vc5'; '0196crit_vc5_md5'];

  age_vs_vel_plots(gauss_fields(idx_4_0G,:), rot_vels, true, 2, [3.0e9, 1.0e10], 'northwest', 'Rotational velocity - 4.0G & var. rotational velocity', 'rot_vel_var_vel_4_0g_z1', num2str(idx));
end


function plot_vel_rot_4_5G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_4_5G;

  age_vs_vel_plots(gauss_fields(idx_4_5G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northwest', 'Rotational velocity - 4.5G & var. rotational velocity', 'rot_vel_var_vel_4_5g', num2str(idx));
end

function plot_vel_rot_4_5G_var_vel2(rot_vels, idx)
  global gauss_fields;
  global idx_4_5G;

  age_vs_vel_plots(gauss_fields(idx_4_5G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northwest', 'Rotational velocity - 4.5G & var. rotational velocity', 'rot_vel_var_vel_4_5g', num2str(idx));
end


function plot_vel_rot_5G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_5_0G;

  age_vs_vel_plots(gauss_fields(idx_5_0G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northwest', 'Rotational velocity - 5.0G & var. rotational velocity', 'rot_vel_var_vel_5_0g', num2str(idx));
end

function plot_vel_rot_5_5G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_5_5G;

  age_vs_vel_plots(gauss_fields(idx_5_5G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northwest', 'Rotational velocity - 5.5G & var. rotational velocity', 'rot_vel_var_vel_5_5g', num2str(idx));
end

function plot_vel_rot_XG_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_X_G;

  age_vs_vel_plots(gauss_fields(idx_X_G,:), rot_vels, true, 20, [1.0e5,1.0e10], 'northwest', '', 'rot_vel_var_vel_var_g', num2str(idx));
end

function plot_vel_rot_XG_var_vel_z1(rot_vels, idx)
  global gauss_fields;
  global idx_X_G;

  age_vs_vel_plots(gauss_fields(idx_X_G,:), rot_vels, true, 2, [3.0e9, 1.0e10], 'northeast', '', 'rot_vel_var_vel_var_g_z1', num2str(idx));
end



function plot_vel_rot_0084vc_var_g(mag_fields, idx)
  global rotational_vels;
  global idx_0084crit;

  age_vs_vel_plots(mag_fields, rotational_vels(idx_0084crit,:), false, 10, [1.0e5,1.0e10], 'northwest', 'Rotational velocity - vcrit=0.0084 & var. magnetic field', 'rot_vel_vc_0084_var_g', num2str(idx));
end

function plot_omega_0_0G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_0_0G;

  age_vs_omega_plots(gauss_fields(idx_0_0G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northeastoutside', 'Omega - 0.0G & var. rotational velocity', 'omega_var_vel_0_0g_', num2str(idx));
end


function plot_omega_3_0G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_3_0G;

  age_vs_omega_plots(gauss_fields(idx_3_0G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northeastoutside', 'Omega - 3.0G & var. rotational velocity', 'omega_var_vel_3_0g', num2str(idx));
end

function plot_omega_4_0G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_4_0G;

  age_vs_omega_plots(gauss_fields(idx_4_0G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northeastoutside', 'Omega - 4.0G & var. rotational velocity', 'omega_var_vel_4_0g', num2str(idx));
end

function plot_omega_XG_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_X_G;

  age_vs_omega_plots(gauss_fields(idx_X_G,:), rot_vels, true, 10, [1.0e5,1.0e10], 'northeastoutside', '', 'omega_var_vel_var_g', num2str(idx));
end



function plot_hr_XG_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_X_G;

  hr_plots(gauss_fields(idx_X_G,:), rot_vels, true, [0.05,0.5], [3.58, 3.8, -0.5, 2.2], 'northwest', '', 'hr_var_vel_var_g', num2str(idx));
end

function plot_hr_XG_var_vel_z1(rot_vels, limits, idx)
  global gauss_fields;
  global idx_X_G;

  hr_plots(gauss_fields(idx_X_G,:), rot_vels, true, [0.005,0.1], limits, 'southwest', '', 'hr_var_vel_var_g_z1', num2str(idx));
end


function plot_hr_0G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_0_0G;

  hr_plots(gauss_fields(idx_0_0G,:), rot_vels, true, [0.05,0.5], [3.58, 3.8, -0.5, 2.2], 'northwest', '', 'hr_var_vel_0_0g_', num2str(idx));
end

function plot_hr_0G_var_vel_z1(rot_vels, idx)
  global gauss_fields;
  global idx_0_0G;

  hr_plots(gauss_fields(idx_0_0G,:), rot_vels, true, [0.005,0.1], [3.75, 3.775, -0.3, 0.45], 'southwest', '', 'hr_var_vel_0_0g_z1', num2str(idx));
end

function plot_hr_2_5G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_2_5G;

  hr_plots(gauss_fields(idx_2_5G,:), rot_vels, true, [0.05,0.5], [3.58, 3.8, -0.5, 2.2], 'southwest', 'HR - 2.5G & var. rotational velocity', 'hr_var_vel_2_5g', num2str(idx));
end


function plot_hr_3G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_3_0G;

  hr_plots(gauss_fields(idx_3_0G,:), rot_vels, true, [0.05,0.5], [3.58, 3.8, -0.5, 2.2], 'southwest', 'HR - 3.0G & var. rotational velocity', 'hr_var_vel_3_0g', num2str(idx));
end


function plot_hr_3_5G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_3_5G;

  hr_plots(gauss_fields(idx_3_5G,:), rot_vels, true, [0.05,0.5], [3.58, 3.8, -0.5, 2.2], 'southwest', 'HR - 3.5G & var. rotational velocity', 'hr_var_vel_3_5g', num2str(idx));
end

function plot_hr_3_5G_var_vel_z_1(rot_vels, idx)
  global gauss_fields;
  global idx_3_5G;

  hr_plots(gauss_fields(idx_3_5G,:), rot_vels, true, [0.005,0.1], [3.75, 3.775, -0.25, 0.4], 'southwest', 'HR - 3.5G & var. rotational velocity', 'hr_var_vel_3_5g_z1', num2str(idx));
end

function plot_hr_4_5G_var_vel2(rot_vels, idx)
  global gauss_fields;
  global idx_4_5G;

  hr_plots(gauss_fields(idx_4_5G,:), rot_vels, true, [0.05,0.5], [3.58, 3.8, -0.5, 2.2], 'southwest', 'HR - 4.5G & var. rotational velocity', 'hr_var_vel_5_0g', num2str(idx));
end


function plot_hr_5_0G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_5_0G;

  hr_plots(gauss_fields(idx_5_0G,:), rot_vels, true, [0.05,0.5], [3.58, 3.8, -0.5, 2.2], 'southwest', 'HR - 5.0G & var. rotational velocity', 'hr_var_vel_5_0g', num2str(idx));
end

function plot_hr_0336vc_var_g(mag_fields, idx)
  global rotational_vels;
  global idx_0336crit;

  hr_plots(mag_fields, rotational_vels(idx_0336crit,:), false, [0.05,0.5], [3.58, 3.8, -0.5, 2.2], 'southwest', 'HR - vcrit=0.0336 & var. magnetic field', 'hr_vc_0336_var_g', num2str(idx));
end

function plot_hr_0336vc_var_g_z1(mag_fields, idx)
  global rotational_vels;
  global idx_0336crit;

  hr_plots(mag_fields, rotational_vels(idx_0336crit,:), false, [0.005,0.1], [3.745, 3.775, -0.3, 0.45], 'southwest', 'HR - vcrit=0.0336 & var. magnetic field', 'hr_vc_0336_var_g_z1', num2str(idx));
end


function plot_cz_size_XG_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_X_G;

  age_vs_cz_size_plots(gauss_fields(idx_X_G,:), rot_vels, true, 0.1, [1.0e5, 1.0e10, 0.0, 1.05], 'northeast', '', 'cz_var_vel_var_g', num2str(idx));
end

function plot_cz_size_XG_var_vel_z1(rot_vels, idx)
  global gauss_fields;
  global idx_X_G;
  age_vs_cz_size_plots(gauss_fields(idx_X_G,:), rot_vels, true, 0.01, [1.0e7, 1.0e10, 0.25, 0.40], 'northeastoutside', '', 'cz_var_vel_var_g_z1', num2str(idx));
end



function plot_cz_size_0G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_0_0G;

  age_vs_cz_size_plots(gauss_fields(idx_0_0G,:), rot_vels, true, 0.1, [1.0e5, 1.0e10, 0.0, 1.05], 'northeast', 'Convective zone radius - 0G & var. rotational velocity', 'cz_var_vel_0g_', num2str(idx));
end

function plot_cz_size_2_5G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_2_5G;

  age_vs_cz_size_plots(gauss_fields(idx_2_5G,:), rot_vels, true, 0.1, [1.0e2, 1.0e10, 0.0, 1.05], 'southwest', 'Convective zone size - 2.5G & var. rotational velocity', 'cz_var_vel_2_5g', num2str(idx));
end


function plot_cz_size_3G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_3_0G;

  age_vs_cz_size_plots(gauss_fields(idx_3_0G,:), rot_vels, true, 0.1, [1.0e2, 1.0e10, 0.0, 1.05], 'southwest', 'Convective zone size - 3.0G & var. rotational velocity', 'cz_var_vel_3_0g', num2str(idx));
end


function plot_cz_size_3_5G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_3_5G;

  age_vs_cz_size_plots(gauss_fields(idx_3_5G,:), rot_vels, true, 0.1, [1.0e2, 1.0e10, 0.0, 1.05], 'southwest', 'Convective zone size - 3.5G & var. rotational velocity', 'cz_var_vel_3_5g', num2str(idx));
end

function plot_cz_size_4_0G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_4_0G;

  age_vs_cz_size_plots(gauss_fields(idx_4_0G,:), rot_vels, true, 0.1, [1.0e2, 1.0e10, 0.0, 1.05], 'southwest', 'Convective zone size - 4.0G & var. rotational velocity', 'cz_var_vel_4_0g', num2str(idx));
end

function plot_cz_size_4_5G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_4_5G;

  age_vs_cz_size_plots(gauss_fields(idx_4_5G,:), rot_vels, true, 0.1, [1.0e2, 1.0e10, 0.0, 1.05], 'southwest', 'Convective zone size - 4.5G & var. rotational velocity', 'cz_var_vel_4_5g', num2str(idx));
end

function plot_cz_size_4_5G_var_vel2(rot_vels, idx)
  global gauss_fields;
  global idx_4_5G;

  age_vs_cz_size_plots(gauss_fields(idx_4_5G,:), rot_vels, true, 0.1, [1.0e2, 1.0e10, 0.0, 1.05], 'southwest', 'Convective zone size - 4.5G & var. rotational velocity', 'cz_var_vel_4_5g', num2str(idx));
end

function plot_cz_size_5_0G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_5_0G;

  age_vs_cz_size_plots(gauss_fields(idx_5_0G,:), rot_vels, true, 0.1, [1.0e2, 1.0e10, 0.0, 1.05], 'southwest', 'Convective zone size - 5.0G & var. rotational velocity', 'cz_var_vel_5_0g', num2str(idx));
end

function plot_cz_size_5_5G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_5_5G;

  age_vs_cz_size_plots(gauss_fields(idx_5_5G,:), rot_vels, true, 0.1, [1.0e2, 1.0e10, 0.0, 1.05], 'southwest', 'Convective zone size - 5.5G & var. rotational velocity', 'cz_var_vel_5_5g', num2str(idx));
end

function plot_cz_size_028vc_var_g(mag_fields, idx)
  global rotational_vels;
  global idx_028crit;

  age_vs_cz_size_plots(mag_fields, rotational_vels(idx_028crit,:), false, 0.1, [1.0e2, 1.0e10, 0.0, 1.05], 'southwest', 'Convective zone radius - vcrit=0.028 & var. magnetic field', 'cz_vc_028_var_g', num2str(idx));
end

function plot_cz_size_028vc_var_g_z1(mag_fields, idx)
  global rotational_vels;
  global idx_028crit;

  age_vs_cz_size_plots(mag_fields, rotational_vels(idx_028crit,:), false, 0.01, [1.0e7, 1.0e10, 0.25, 0.30], 'north', 'Convective zone radius - vcrit=0.028 & var. magnetic field','cz_vc_028_var_g_z1', num2str(idx));
end

function plot_cz_size_028vc_var_g_z1_special(mag_fields, idx)
  global rotational_vels;
  global idx_028crit;

  age_vs_cz_size_plots_special(rotational_vels(idx_028crit,:), false, 0.01, [1.0e7, 1.0e10, 0.25, 0.30], 'north', 'Convective zone radius - vcrit=0.028 & var. magnetic field','cz_vc_028_var_g_z1', num2str(idx));
end

function plot_cz_size_0G_var_vel_z1(rot_vels,idx)
  global gauss_fields;
  global idx_0_0G;

  age_vs_cz_size_plots(gauss_fields(idx_0_0G,:), rot_vels, true, 0.01, [1.0e7, 1.0e10, 0.25, 0.30], 'north', 'Convective zone radius - 0G & var. rotational velocity', 'cz_var_vel_0_0g_z1', num2str(idx));
end



function plot_cz_size_3_5G_var_vel_z1(rot_vels, idx)
  global gauss_fields;
  global idx_3_5G;

  age_vs_cz_size_plots(gauss_fields(idx_3_5G,:), rot_vels, true, 0.01, [1.0e7, 1.0e10, 0.25, 0.28], 'north', 'Convective zone size - 3.5G & var. rotational velocity', 'cz_var_vel_3_5g_z1', num2str(idx));
end


function plot_m_dot_028vc_var_g(mag_fields, idx)
  global rotational_vels;
  global idx_028crit;

  age_vs_m_dot_plots(mag_fields, rotational_vels(idx_028crit,:), false, 0.5, [1.0e2, 1.0e10, -14.0, -10.0], 'southwest', 'Mass loss - vcrit=0.028 & var. magnetic field', 'mdot_vc_028_var_g', num2str(idx));
end

function plot_m_dot_028vc_var_g_z1(mag_fields, idx)
  global rotational_vels;
  global idx_028crit;

  age_vs_m_dot_plots(mag_fields, rotational_vels(idx_028crit,:), false, 0.05, [3.0e8, 1.0e10, -13.6, -13.3], 'northwest', 'Mass loss - vcrit=0.028 & var. magnetic field', 'mdot_vc_028_var_g_z1', num2str(idx));
end

function plot_m_dot_0G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_0_0G;

  age_vs_m_dot_plots(gauss_fields(idx_0_0G,:), rot_vels, true, 0.5, [1.0e2, 1.0e10, -14.0, -10.0], 'southwest', 'Mass loss - 0.0G & var. rotational velocity', 'mdot_var_vel_0_0_g_', num2str(idx));
end

function plot_m_dot_3G_var_vel(rot_vels, idx)
  global gauss_fields;
  global idx_3_0G;

  age_vs_m_dot_plots(gauss_fields(idx_3_0G,:), rot_vels, true, 0.5, [1.0e2, 1.0e10, -14.0, -10.0], 'southwest', 'Mass loss - 3.0G & var. rotational velocity', 'mdot_var_vel_3_0_g', num2str(idx));
end

function plot_m_dot_XG_var_vel(rot_vels,idx)
  global gauss_fields;
  global idx_X_G;

  age_vs_m_dot_plots(gauss_fields(idx_X_G,:), rot_vels, true, 0.5, [1.0e5, 1.0e10, -14.0, -10.0], 'northwest', '', 'mdot_var_vel_g', num2str(idx));
end

function plot_m_dot_XG_var_vel_z1(rot_vels,idx)
  global gauss_fields;
  global idx_X_G;

  age_vs_m_dot_plots(gauss_fields(idx_X_G,:), rot_vels, true, 0.05, [3.0e8, 1.0e10, -13.6, -13.3], 'northwest', '', 'mdot_var_vel_g_z1', num2str(idx));
end



function plot_reimers_m_dot_0G_var_vel(rot_vels)
  global gauss_fields;
  global idx_0_0G;

  age_vs_reimers_m_dot_plots(gauss_fields(idx_0_0G,:), rot_vels, true, 0.5, [1.0e2, 1.0e10, -0.0, 2.5], 'Reimers mass loss - 0.0G & var. rotational velocity', 'mdot_var_vel_0_0_g');
end

function plot_m_l_r_0G_var_vel()
  global gauss_fields;
  global idx_0_0G;
  global rotational_vels;

  age_vs_m_l_r_plots(gauss_fields(idx_0_0G,:), rotational_vels(1:1,:), true, 5.0, [1.0e2, 1.0e10, 0.0, 100], 'M, L, R evolution - 0.0G & var. rotational velocity', 'mdot_var_vel_0_0_g');
end

function plot_m_l_r_0G_var_vel_z1()
  global gauss_fields;
  global idx_0_0G;
  global rotational_vels;

  age_vs_m_l_r_plots(gauss_fields(idx_0_0G,:), rotational_vels(1:1,:), true, 0.1, [1.0e6, 1.0e10, 0.0, 100], 'M, L, R evolution - 0.0G & var. rotational velocity', 'mdot_var_vel_0_0_g');
end



function plot_age_vs_mb_activation_2_5G(rot_vels, idx)
  global gauss_fields;
  global idx_2_5G;

  age_vs_mb_activation(gauss_fields(idx_2_5G,:), rot_vels, true, 'eastoutside', 'Magnetic braking activation & core nature - 2.5G & var. rotational velocity', 'mb_act_var_vel_2_5g', num2str(idx));
end


function plot_age_vs_mb_activation_3_0G(rot_vels, idx)
  global gauss_fields;
  global idx_3_0G;

  age_vs_mb_activation(gauss_fields(idx_3_0G,:), rot_vels, true, 'eastoutside', 'Magnetic braking activation & core nature - 3.0G & var. rotational velocity', 'mb_act_var_vel_3_0g', num2str(idx));
end


function plot_age_vs_mb_activation_3_5G(rot_vels, idx)
  global gauss_fields;
  global idx_3_5G;

  age_vs_mb_activation(gauss_fields(idx_3_5G,:), rot_vels, true, 'eastoutside', 'Magnetic braking activation & core nature - 3.5G & var. rotational velocity', 'mb_act_var_vel_3_5g', num2str(idx));
end

function plot_age_vs_mb_activation_4_0G(rot_vels, idx)
  global gauss_fields;
  global idx_4_0G;

  age_vs_mb_activation(gauss_fields(idx_4_0G,:), rot_vels, true, 'eastoutside', 'Magnetic braking activation & core nature - 4.0G & var. rotational velocity', 'mb_act_var_vel_4_0g', num2str(idx));
end

function plot_age_vs_mb_activation_4_5G(rot_vels, idx)
  global gauss_fields;
  global idx_4_5G;

  age_vs_mb_activation(gauss_fields(idx_4_5G,:), rot_vels, true, 'eastoutside', 'Magnetic braking activation & core nature - 4.5G & var. rotational velocity', 'mb_act_var_vel_4_5g', num2str(idx));
end

function plot_age_vs_mb_activation_5_0G(rot_vels, idx)
  global gauss_fields;
  global idx_5_0G;

  age_vs_mb_activation(gauss_fields(idx_5_0G,:), rot_vels, true, 'eastoutside', 'Magnetic braking activation & core nature - 5.0G & var. rotational velocity', 'mb_act_var_vel_5_0g', num2str(idx));
end

function plot_age_vs_mb_activation_5_5G(rot_vels, idx)
  global gauss_fields;
  global idx_5_5G;

  age_vs_mb_activation(gauss_fields(idx_5_5G,:), rot_vels, true, 'eastoutside', 'Magnetic braking activation & core nature - 5.5G & var. rotational velocity', 'mb_act_var_vel_5_5g', num2str(idx));
end

function plot_age_vs_mb_activation_XG(rot_vels,idx)
  global gauss_fields;
  global idx_X_G;

  age_vs_mb_activation(gauss_fields(idx_X_G,:), rot_vels, true, 'northeastoutside', '', 'mb_act_var_vel_g', num2str(idx));
end


function plot_age_vs_mb_activation_0084vc(mag_fields, idx)
  global rotational_vels;
  global idx_0084crit;

  age_vs_mb_activation(mag_fields, rotational_vels(idx_0084crit,:), false, 'eastoutside', 'Magnetic braking activation & core nature - vcrit=0.0084 & var. magnetic field', 'mb_act_vc_0084_var_g', num2str(idx));
end

function plot_age_vs_mb_activation_028vc(mag_fields, idx)
  global rotational_vels;
  global idx_028crit;

  age_vs_mb_activation(mag_fields, rotational_vels(idx_028crit,:), false, 'eastoutside', 'MB activation & core nature - vcrit=0.028 & var. magnetic field', 'mb_act_vc_028_var_g', num2str(idx));
end


function plot_age_vs_alpha_mlt_3_0G(rot_vels, idx)
  global gauss_fields;
  global idx_3_0G;

  age_vs_alpha_mlt_plots(gauss_fields(idx_3_0G,:), rot_vels, true, 0.02, [1.0e2, 1.0e10, 1.65, 1.9], 'eastoutside', '\alpha_{MLT} - 3.0G & var. rotational velocity', 'alpha_mlt_var_vel_3_0g', num2str(idx));
end

function plot_age_vs_alpha_mlt_4_0G(rot_vels, idx)
  global gauss_fields;
  global idx_4_0G;

  age_vs_alpha_mlt_plots(gauss_fields(idx_4_0G,:), rot_vels, true, 0.02, [1.0e2, 1.0e10, 1.65, 1.9], 'eastoutside', '\alpha_{MLT} - 4.0G & var. rotational velocity', 'alpha_mlt_var_vel_4_0g', num2str(idx));
end

function plot_age_vs_alpha_mlt_XG(rot_vels, idx)
  global gauss_fields;
  global idx_X_G;

  age_vs_alpha_mlt(gauss_fields(idx_X_G,:), rot_vels, true, 0.02, [1.0e5, 1.0e10, 1.65, 1.88], 'northwest', '', 'alpha_mlt_var_vel_g', num2str(idx));
end

function plot_omega_vs_mag_field_XG(rot_vels, show_limits, idx)
  global gauss_fields;
  global idx_X_G;

  %omegs_vs_mag_field_plots(gauss_fields(idx_X_G,:), rot_vels, true, 10, [0.5, 50, 1, 1000], 'northwest', 'Magnetic field intensity  & var. rotational velocity', 'mag_field_var_vel_g');
  omegs_vs_mag_field(gauss_fields(idx_X_G,:), rot_vels, show_limits, 10, [1.0e5, 9.0e9, 1, 3000], 'southwest', '', 'mag_field_var_vel_g', num2str(idx));
end

function plot_age_vs_teff_XG(rot_vels, idx)
  global gauss_fields;
  global idx_X_G;

  age_vs_teff_plots(gauss_fields(idx_X_G,:), rot_vels, true, 0.5, [1.0e5, 1.0e10, 2.5, 5],  'northwest', '', 'teff_logg_var_vel_g', num2str(idx));
end

function plot_age_vs_teff_XG_z1(rot_vels, idx)
  global gauss_fields;
  global idx_X_G;

  age_vs_teff_plots(gauss_fields(idx_X_G,:), rot_vels, true, 0.2, [1.0e6, 1.0e9, 3.6, 4.6],  'northwest', '', 'teff_logg_var_vel_g_z1', num2str(idx));
end

function plot_radius_vs_mag_field_XG(rot_vels, idx)
  global gauss_fields;
  global idx_X_G;

  age_vs_radius_plots(gauss_fields(idx_X_G,:), rot_vels, true, 0.2, [1.0e5, 1.0e10, -0.2, 0.6],  'northwest', '', 'lograd_var_vel_g', num2str(idx));
end

function plot_radius_vs_mag_field_XG_z1(rot_vels, idx)
  global gauss_fields;
  global idx_X_G;

  age_vs_radius_plots(gauss_fields(idx_X_G,:), rot_vels, true, 0.05, [1.0e6, 1.0e9, -0.1, 0.1],  'northwest', '', 'lograd_var_vel_g_z1', num2str(idx));
end



function plot_kipperhahn_3G_var_vel(rot_vels)
  global gauss_fields;
  global idx_3_0G;

  kipperhahn_plots(gauss_fields(idx_3_0G,:), rot_vels, true, 0.1, [1.0e2, 1.0e10, 0.0, 1.05], 'southwest', '', 'kp_var_vel_3_0g');
end


function paper1()
  global gauss_fields;
  global rotational_vels;

  global idx_0_0G;
  global idx_3_0G;
  global idx_3_5G;
  global idx_4_0G;
  global idx_4_5G;
  global idx_5_0G;

  global idx_0crit;
  global idx_0084crit;
  global idx_014crit;
  global idx_0196crit;
  global idx_028crit;
  global idx_0336crit;
  global idx_029crit;
  global idx_030crit;
  global idx_031crit;
  global idx_032crit;
  global idx_0314crit;


  mag_fields = gauss_fields([idx_0_0G;idx_3_0G;idx_3_5G;idx_4_0G;idx_4_5G;idx_5_0G],:);
  mag_fields2 = gauss_fields([idx_3_0G;idx_3_5G;idx_4_0G;idx_4_5G;idx_5_0G],:);

  rot_vels = rotational_vels([idx_0crit;idx_0084crit;idx_014crit;idx_0196crit;idx_028crit;idx_0336crit],:);
  rot_vels2 = rotational_vels([idx_0084crit;idx_014crit;idx_0196crit;idx_028crit;idx_0336crit],:);
  rot_vels3 = rotational_vels([idx_0314crit],:);


  plot_0G_var_vel(rot_vels, 0);
  plot_0G_var_vel_z1(rot_vels, 0);
  plot_vel_rot_0G_var_vel(rot_vels2, 0);
  %plot_hr_0336vc_var_g();
  plot_hr_0336vc_var_g_z1(mag_fields, 0336);
  plot_hr_0G_var_vel_z1(rot_vels, 0);
  plot_4_0G_var_vel(rot_vels2, 4);
  plot_vel_rot_4G_var_vel(rot_vels2, 4);
  plot_vel_rot_4G_var_vel_z1(rot_vels2, 4);
  plot_cz_size_028vc_var_g(mag_fields, 028);
  plot_cz_size_028vc_var_g_z1(mag_fields, 028);
  plot_m_dot_028vc_var_g(mag_fields, 028);
  plot_m_dot_028vc_var_g_z1(mag_fields, 028);
  plot_age_vs_mb_activation_028vc(mag_fields2, 028);
  plot_3_0G_0314vc(rot_vels3, 3);
  plot_vel_rot_3_0G_0314vc(rot_vels3, 3);

  %grid Li var_vel
  plot_0G_var_vel(rot_vels, 0);
  plot_3_0G_var_vel(rot_vels2, 3);
  plot_3_5G_var_vel(rot_vels2, 35);
  plot_4_0G_var_vel(rot_vels2, 4);
  plot_4_5G_var_vel(rot_vels2, 45);
  plot_5_0G_var_vel(rot_vels2, 5);
  %plot_5_5G_var_vel();

  %grid Li var_b
  plot_0084vc_var_g(mag_fields, 0084);
  plot_014vc_var_g(mag_fields, 014);
  plot_0196vc_var_g(mag_fields, 0196);
  plot_028vc_var_g(mag_fields, 028);
  plot_0336vc_var_g(mag_fields, 0336);

  %grid rot vel
  plot_vel_rot_0G_var_vel(rot_vels2, 0);
  plot_vel_rot_3_0G_var_vel(rot_vels2, 3);
  plot_vel_rot_3_5G_var_vel(rot_vels2, 35);
  plot_vel_rot_4G_var_vel(rot_vels2, 4);
  plot_vel_rot_4_5G_var_vel(rot_vels2, 45);
  plot_vel_rot_5G_var_vel(rot_vels2, 5);
  %plot_vel_rot_5_5G_var_vel();

  %grid mag breaking
  %plot_age_vs_mb_activation_3_0G(rot_vels2);
  %plot_age_vs_mb_activation_3_5G(rot_vels2);
  %plot_age_vs_mb_activation_4_0G(rot_vels2);
  %plot_age_vs_mb_activation_4_5G(rot_vels2);
  %plot_age_vs_mb_activation_5_0G(rot_vels2);
  %plot_age_vs_mb_activation_5_5G(rot_vels2);

end



function paper2()
  global gauss_fields;
  global rotational_vels;
  global idx_0crit;
  global idx_0084crit;
  global idx_014crit;
  global idx_0196crit;
  global idx_028crit;
  global idx_0336crit;
  global idx_11crit;
  global idx_105crit;
  global idx_1075crit;
  global idx_1125crit;
  global idx_1025crit;
  global idx_115crit;
  global idx_10crit;
  global idx_0975crit;
  global idx_095crit;
  global idx_0925crit;
  global idx_125crit;
  global idx_135crit;
  global idx_1175crit;
  global idx_12crit;
  global idx_1225crit;
  global idx_1275crit;
  global idx_13crit;
  global idx_1325crit;
  global idx_135crit;
  global idx_1375crit;
  global idx_14crit;
  global idx_1425crit;
  global idx_145crit;
  global idx_1475crit;
  global idx_15crit;
  global idx_1525crit;
  global idx_155crit;


  #rot_vels5 = rotational_vels([idx_0925crit;idx_095crit;idx_0975crit;idx_10crit;idx_1025crit],:);
  #rot_vels6 = rotational_vels([idx_105crit;idx_1075crit;idx_11crit;idx_1125crit;idx_115crit],:);
  #rot_vels7 = rotational_vels([idx_1175crit;idx_12crit;idx_1225crit;idx_125crit;idx_1275crit],:);
  #rot_vels8 = rotational_vels([idx_13crit;idx_1325crit;idx_1375crit;idx_14crit;idx_1425crit],:);

  rot_vels = rotational_vels([idx_0crit;idx_0084crit;idx_014crit;idx_0196crit;idx_028crit;idx_0336crit],:);
  rot_vels5 = rotational_vels([idx_095crit;idx_10crit;idx_105crit;idx_11crit;idx_115crit],:);
  #rot_vels5 = rotational_vels([idx_095crit],:);
  rot_vels6 = rotational_vels([idx_11crit;idx_115crit;idx_12crit],:);
  rot_vels7 = rotational_vels([idx_12crit;idx_125crit;idx_13crit;idx_14crit;idx_1425crit;],:);
  rot_vels8 = rotational_vels([idx_13crit;idx_14crit;idx_1425crit],:);


  #plot_XG_var_vel(rotational_vels([idx_1425crit],:),1);
  plot_0G_var_vel(rot_vels,0);
  plot_4_0G_var_vel(rot_vels,4);
  #plot_0G_var_vel_z1(rot_vels,0);
  plot_XG_var_vel(rot_vels5,1);
  #plot_XG_var_vel(rot_vels6,2);
  plot_XG_var_vel(rot_vels7,3);
  #plot_XG_var_vel(rot_vels8,4);

  plot_age_vs_alpha_mlt_XG(rot_vels5,1);
  #plot_age_vs_alpha_mlt_XG(rot_vels6,2);
  plot_age_vs_alpha_mlt_XG(rot_vels7,3);
  #plot_age_vs_alpha_mlt_XG(rot_vels8,4);

  plot_vel_rot_0G_var_vel(rot_vels,0);
  plot_vel_rot_XG_var_vel(rot_vels5,1);
  #plot_vel_rot_XG_var_vel(rot_vels6,2);
  plot_vel_rot_XG_var_vel(rot_vels7,3);
  #plot_vel_rot_XG_var_vel(rot_vels8,4);

  #plot_omega_XG_var_vel(rot_vels8,3);

  plot_omega_vs_mag_field_XG(rot_vels5, false, 1);
  #plot_omega_vs_mag_field_XG(rot_vels6, false, 2);
  plot_omega_vs_mag_field_XG(rot_vels7, false, 3);
  #plot_omega_vs_mag_field_XG(rot_vels8, false, 4);
  #plot_omega_vs_mag_field_XG(rotational_vels([idx_1425crit],:), true, 3);

  #plot_cz_size_0G_var_vel_z1(rot_vels,0);
  plot_cz_size_XG_var_vel_z1(rot_vels5,1);
  #plot_cz_size_XG_var_vel_z1(rot_vels6,2);
  plot_cz_size_XG_var_vel_z1(rot_vels7,3);
  #plot_cz_size_XG_var_vel_z1(rot_vels8,4);

  plot_cz_size_0G_var_vel(rot_vels,0);
  plot_cz_size_XG_var_vel(rot_vels5,1);
  #plot_cz_size_XG_var_vel(rot_vels6,2);
  plot_cz_size_XG_var_vel(rot_vels7,3);
  #plot_cz_size_XG_var_vel(rot_vels8,4);

  plot_hr_0G_var_vel(rot_vels,0);
  plot_hr_XG_var_vel(rot_vels5,1);
  #plot_hr_XG_var_vel(rot_vels6,2);
  plot_hr_XG_var_vel(rot_vels7,3);
  #plot_hr_XG_var_vel(rot_vels8,4);

  plot_hr_0G_var_vel_z1(rot_vels,0);
  plot_hr_XG_var_vel_z1(rot_vels5,[3.70, 3.73, -0.35, 0.0],1);
  #plot_hr_XG_var_vel_z1(rot_vels6,[3.70, 3.73, -0.35, 0.0],2);
  plot_hr_XG_var_vel_z1(rot_vels7,[3.68, 3.72, -0.45, -0.10],3);
  #plot_hr_XG_var_vel_z1(rot_vels8,[3.68, 3.71, -0.45, -0.15],4);

  #plot_age_vs_mb_activation_XG(rot_vels5,1);
  #plot_age_vs_mb_activation_XG(rot_vels6,2);
  plot_age_vs_mb_activation_XG(rot_vels7,3);
  #plot_age_vs_mb_activation_XG(rot_vels8,4);

  plot_m_dot_XG_var_vel(rot_vels5,1);
  #plot_m_dot_XG_var_vel(rot_vels6,2);
  plot_m_dot_XG_var_vel(rot_vels7,3);
  #plot_m_dot_XG_var_vel(rot_vels8,4);

  plot_m_dot_XG_var_vel_z1(rot_vels5,1);
  #plot_m_dot_XG_var_vel_z1(rot_vels6,2);
  plot_m_dot_XG_var_vel_z1(rot_vels7,3);
  #plot_m_dot_XG_var_vel_z1(rot_vels8,4);

  plot_age_vs_teff_XG(rot_vels5,1);
  #plot_age_vs_teff_XG(rot_vels6,2);
  plot_age_vs_teff_XG(rot_vels7,3);
  #plot_age_vs_teff_XG(rot_vels8,4);

  plot_age_vs_teff_XG_z1(rot_vels5,1);
  #plot_age_vs_teff_XG(rot_vels6,2);
  plot_age_vs_teff_XG_z1(rot_vels7,3);
  #plot_age_vs_teff_XG(rot_vels8,4);

  plot_radius_vs_mag_field_XG(rot_vels5,1);
  #plot_radius_vs_mag_field_XG(rot_vels6,2);
  plot_radius_vs_mag_field_XG(rot_vels7,3);
  #plot_radius_vs_mag_field_XG(rot_vels8,4);

  plot_radius_vs_mag_field_XG_z1(rot_vels5,1);
  #plot_radius_vs_mag_field_XG(rot_vels6,2);
  plot_radius_vs_mag_field_XG_z1(rot_vels7,3);
  #plot_radius_vs_mag_field_XG(rot_vels8,4);

end


function main()
  global gauss_fields;
  global rotational_vels;
  global dl_rotational_vels;

  global idx_0_0G;
  global idx_3_0G;
  global idx_3_5G;
  global idx_4_0G;
  global idx_4_5G;
  global idx_5_0G;

  global idx_0crit;
  global idx_0084crit;
  global idx_014crit;
  global idx_0196crit;
  global idx_028crit;
  global idx_0336crit;
  global idx_029crit;
  global idx_030crit;
  global idx_031crit;
  global idx_032crit;
  global idx_0312crit;
  global idx_0314crit;
  global idx_9_090256e_6_dl;
  global idx_0336crit_alpha;
  global idx_11crit;
  global idx_105crit;
  global idx_1075crit;
  global idx_1125crit;
  global idx_1025crit;
  global idx_115crit;
  global idx_10crit;
  global idx_0975crit;
  global idx_095crit;
  global idx_0925crit;
  global idx_125crit;
  global idx_135crit;
  global idx_1175crit;
  global idx_12crit;
  global idx_1225crit;
  global idx_1275crit;
  global idx_13crit;
  global idx_1325crit;
  global idx_1375crit;
  global idx_14crit;
  global idx_1425crit;
  global idx_145crit;
  global idx_1475crit;
  global idx_15crit;
  global idx_1525crit;
  global idx_155crit;



  mag_fields = gauss_fields([idx_0_0G;idx_3_0G;idx_3_5G;idx_4_0G;idx_4_5G;idx_5_0G],:);
  mag_fields2 = gauss_fields([idx_3_0G;idx_3_5G;idx_4_0G;idx_4_5G;idx_5_0G],:);
  mag_fields3 = gauss_fields([idx_3_0G],:);
  mag_fields_annexA = gauss_fields([idx_5_0G],:);
  mag_fields_annexB = gauss_fields([idx_5_0G],:);

  rot_vels = rotational_vels([idx_0084crit;idx_014crit;idx_0196crit;idx_028crit;idx_0336crit],:);
  dl_rot_vels = dl_rotational_vels([idx_0336crit],:);
  rot_vels2 = rotational_vels([idx_0084crit;idx_014crit;idx_0196crit;idx_028crit;idx_0336crit],:);
  %rot_vels3 = rotational_vels([idx_030crit;idx_031crit;idx_0312crit;idx_0314crit;idx_032crit],:);
  rot_vels3 = rotational_vels([idx_0314crit],:);
  rot_vels4 = rotational_vels([idx_028crit;idx_0314crit],:);
  rot_vels5 = rotational_vels([idx_0925crit;idx_095crit;idx_0975crit;idx_10crit;idx_1025crit],:);
  rot_vels6 = rotational_vels([idx_105crit;idx_1075crit;idx_11crit;idx_1125crit;idx_115crit],:);
  rot_vels7 = rotational_vels([idx_1175crit;idx_12crit;idx_1225crit;idx_125crit;idx_1275crit],:);
  rot_vels8 = rotational_vels([idx_13crit;idx_1325crit;idx_1375crit;idx_14crit;idx_1425crit],:);
  rot_vels9 = rotational_vels([idx_145crit;idx_1475crit;idx_15crit;idx_1525crit;idx_155crit],:);
  rot_vels10 = rotational_vels([idx_125crit;idx_13crit;idx_14crit;idx_1475crit;idx_155crit],:);
  rot_vels11 = rotational_vels([idx_12crit;idx_125crit;idx_13crit;idx_14crit;idx_1425crit;],:);






  rot_vels_annexB = rotational_vels([idx_028crit],:);

  %plot_age_vs_mb_activation_3_5G();
  %plot_age_vs_mb_activation_4_0G();
  %plot_age_vs_mb_activation_4_5G();
  %plot_age_vs_mb_activation_5_0G();
  %plot_age_vs_mb_activation_5_5G();
  %plot_age_vs_mb_activation_0084vc();
  %plot_age_vs_mb_activation_028vc();
  %plot_age_vs_mb_activation_028vc(mag_fields2);


  %plot_0G_var_vel(rot_vels);
  %plot_0G_var_vel_z1(rot_vels);
  %plot_2_5G_var_vel(rot_vels);
  %plot_3_0G_var_vel(rot_vels);
  %plot_3_0G_var_vel(rot_vels3);
  %plot_3_0G_var_vel(rotational_vels([idx_0336crit],:));

  %plot_3_0G_var_vel(dl_rotational_vels([idx_0336crit],:));

  %plot_3_0G_var_vel(rotational_vels([idx_0336crit],:));  rot_vels6 = rotational_vels([idx_105crit;idx_1075crit;idx_11crit;idx_1125crit;idx_115crit],:);
  %plot_3_0G_var_vel(rotational_vels([idx_0336crit,idx_0336crit_alpha],:));
  %plot_3_0G_0314vc(rot_vels3);
  %plot_3_5G_var_vel(rot_vels);
  %plot_4_0G_var_vel(rot_vels);
  %plot_4_3G_var_vel(rot_vels);
  %plot_4_5G_var_vel(rot_vels);
  %plot_5_0G_var_vel(rot_vels);
  %plot_5_5G_var_vel(rot_vels);

  %Works only with paper3d folder
  %plot_XG_var_vel(rot_vels5);
  %plot_XG_var_vel(rot_vels6);
  %plot_XG_var_vel(rot_vels7);
  %plot_XG_var_vel(rot_vels8);
  %plot_XG_var_vel(rot_vels9);
  %plot_XG_var_vel(rot_vels10);
  %plot_XG_var_vel(rotational_vels([idx_135crit],:));
  %plot_XG_var_vel(rotational_vels([idx_1075crit],:));
  %plot_XG_var_vel(rotational_vels([idx_125crit],:), 1);
  %plot_XG_var_vel(rotational_vels([idx_1475crit],:));

  %plot_age_vs_alpha_mlt_3_0G(dl_rotational_vels([idx_9_090256e_6_dl],:));
  %plot_age_vs_alpha_mlt_3_0G(rotational_vels([idx_0336crit_alpha],:));
  %plot_age_vs_alpha_mlt_4_0G(rot_vels2);
  %plot_age_vs_alpha_mlt_XG(rot_vels6);
  %plot_age_vs_alpha_mlt_XG(rot_vels7);
  %plot_age_vs_alpha_mlt_XG(rotational_vels([idx_1025crit],:));
  %plot_age_vs_alpha_mlt_XG(rotational_vels([idx_1475crit],:),1);
  %plot_age_vs_teff_XG_z1(rot_vels11,1);

  %plot_kipperhahn_3G_var_vel(dl_rotational_vels([idx_0336crit],:));

  %Includes dynamo effect
  %plot_4_0G_var_vel_st();

  %Var control 0.00001
  %plot_4_0G_var_vel_vc5();

  % Pruebas ad hoc
  %plot_4_5G_var_vel_special();
  %plot_vel_rot_4_5G_var_vel2();
  %plot_cz_size_4_5G_var_vel2();
  %plot_hr_4_5G_var_vel2();
  %plot_2G_var_vel();
  %plot_3G_var_vel();
  %plot_3_3G_var_vel();
  %plot_vel_rot_2G_var_vel();
  %plot_vel_rot_3_0G_var_vel(rot_vels3);
  %plot_vel_rot_3_0G_var_vel(rotational_vels([idx_0336crit],:));

  %plot_vel_rot_3_0G_var_vel(dl_rotational_vels([idx_0336crit],:));

  %plot_vel_rot_3_0G_0314vc(rot_vels3);
  %plot_vel_rot_3_0G_var_vel_mlt(rot_vels4);
  %plot_vel_rot_3_5G_var_vel();
  %plot_vel_rot_3_3G_var_vel();

  %plot_omega_0_0G_var_vel(rotational_vels([idx_0336crit],:));
  %plot_omega_3_0G_var_vel(rotational_vels([idx_0336crit],:));
  %plot_omega_3_0G_var_vel(dl_rotational_vels([idx_0336crit],:));
  %plot_omega_0_0G_var_vel(rot_vels2);
  %plot_omega_3_0G_var_vel(rot_vels2);
  %plot_omega_4_0G_var_vel(rot_vels2)

  %plot_0084vc_var_g();
  %plot_014vc_var_g();
  %plot_0196vc_var_g();
  %plot_028vc_var_g();
  %plot_0336vc_var_g();


  %plot_vel_rot_0G_var_vel(rot_vels2);
  %plot_vel_rot_2_5G_var_vel(rot_vels2);
  %plot_vel_rot_3_0G_var_vel(rot_vels2);
  %plot_vel_rot_3_5G_var_vel(rot_vels2);
  %plot_vel_rot_4G_var_vel(rot_vels2);
  %plot_vel_rot_4G_var_vel_z1(rot_vels2);
  %plot_vel_rot_4_5G_var_vel(rot_vels2);
  %plot_vel_rot_5G_var_vel(rot_vels2);
  %plot_vel_rot_5_5G_var_vel(rot_vels2);
  %plot_vel_rot_XG_var_vel(rot_vels2);
  %plot_vel_rot_XG_var_vel_z1(rot_vels2);
  %plot_vel_rot_XG_var_vel(rot_vels6);
  %plot_vel_rot_XG_var_vel(rotational_vels([idx_105crit],:));
  %plot_vel_rot_XG_var_vel(rot_vels7);
  %plot_vel_rot_XG_var_vel(rot_vels8);
  %plot_vel_rot_XG_var_vel(rot_vels9);
  %plot_vel_rot_XG_var_vel(rot_vels10);
  %plot_vel_rot_XG_var_vel(rotational_vels([idx_135crit],:),1);
  %plot_vel_rot_XG_var_vel(rotational_vels([idx_1475crit],:),3);

  %Var control 0.00001
  %plot_vel_rot_4G_var_vel_vc5_md5();
  %plot_vel_rot_4G_var_vel_z1_vc5_md5();


  %plot_cz_size_0G_var_vel();
  %plot_cz_size_3G_var_vel(rot_vels4);
  %plot_cz_size_3G_var_vel(dl_rotational_vels([idx_0336crit],:));
  %plot_cz_size_3_5G_var_vel();
  %plot_cz_size_3_5G_var_vel_z_1();
  %plot_cz_size_4_0G_var_vel();
  %plot_cz_size_4_5G_var_vel();
  %plot_cz_size_5_0G_var_vel();
  %plot_cz_size_5_5G_var_vel();
  %plot_cz_size_028vc_var_g(mag_fields);
  %plot_cz_size_028vc_var_g_z1(mag_fields);
  %plot_cz_size_028vc_var_g_z1_special(mag_fields);
  %plot_cz_size_028vc_var_g_z1_special();set(gca,'YTick',0:ytick:1.0);
  %plot_cz_size_0G_var_vel_z1();
  %plot_cz_size_XG_var_vel(rotational_vels([idx_105crit],:));
  %plot_cz_size_XG_var_vel_z1(rotational_vels([idx_105crit],:));
  %plot_cz_size_XG_var_vel_z1(rot_vels6);


  %plot_m_dot_028vc_var_g();
  %plot_m_dot_028vc_var_g(mag_fields);
  %plot_m_dot_028vc_var_g_z1();
  %plot_m_dot_0G_var_vel();
  %plot_m_dot_3G_var_vel(rot_vels4);
  %plot_m_dot_XG_var_vel(rotational_vels([idx_125crit],:),1)
  %plot_reimers_m_dot_0G_var_vel();
  %plot_m_l_r_0G_var_vel();
  %plot_m_l_r_0G_var_vel_z1();


  %plot_hr_3_0G_var_vel();
  %plot_hr_3G_var_vel(rot_vels4);
  %plot_hr_3G_var_vel(rotational_vels([idx_0336crit],:));
  %plot_hr_3G_var_vel(dl_rotational_vels([idx_0336crit],:));
  %plot_hr_3_5G_var_vel();
  %plot_hr_5_0G_var_vel();
  %plot_hr_0336vc_var_g();
  %plot_hr_0336vc_var_g_z1(mag_fields);
  %plot_hr_0G_var_vel();
  %plot_hr_0G_var_vel_z1();
  %plot_hr_3_5G_var_vel_z_1();
  %plot_hr_XG_var_vel(rot_vels6);
  %plot_hr_XG_var_vel(rot_vels8);
  %plot_hr_XG_var_vel(rotational_vels([idx_135crit],:),1);

  %plot_omega_vs_mag_field_XG(rot_vels6, false);
  %plot_omega_vs_mag_field_XG(rotational_vels([idx_0975crit],:), false);
  %plot_omega_vs_mag_field_XG(rotational_vels([idx_0975crit],:), true);
  %plot_omega_vs_mag_field_XG(rotational_vels([idx_1425crit],:), true, 3);
  %plot_omega_vs_mag_field_XG(rot_vels7, false, 3);
  %plot_0G_var_vel(rot_vels,0);
  %plot_hr_0G_var_vel(rot_vels,0);
  %plot_vel_rot_0G_var_vel(rot_vels,0);
  %plot_hr_0G_var_vel_z1(rot_vels,0);
  %plot_age_vs_teff_XG(rot_vels11,3);
  %plot_age_vs_teff_XG_z1(rot_vels7,3);
  %plot_4_0G_var_vel(rot_vels,4);
  %plot_vel_rot_4G_var_vel(rot_vels,4);
  %plot_vel_rot_4G_var_vel_z1(rot_vels,4);
  %plot_omega_vs_mag_field_XG(rot_vels11, false, 3);

  %filename  = "/home/rcaballeron/MESA/workspace/sun-jupiter-system/Docs/runs/run_paper/4g_12kms/1M_photosphere_history.data";
  %calculate_ZAMS(filename);

  %plot_teff_vs_mag_field_XG(rotational_vels([idx_1475crit],:),1);

  paper1();
  %paper2();
end


main();
