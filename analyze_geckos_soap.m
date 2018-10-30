%% Start up

clc;
clear all;
close all;
path(pathdef);

%% Load all data into a structure array
addpath '/Users/jnirody/Documents/Dropbox/Research/geckos/Data/All trials and outputs/ALL trials/Soap txt files' % where are the files
num_geckos = 10; % highest numbered gecko
num_trials = 23; % highest numbered trial

% measurements (all in mm)

 weights = [0,0,5.4,0,5.22,4.89,5.71,6.9,6.6,5.3];
 sex = ['f','f','f','m','m','m','m','m','m','f'];
 svl = [0,0,59.91,0,0,54.4,59.6,57.4,58.6,56.5];
 l_front = [0,0,6.2,0,6,7,6.5,6.2,7.4,6.6];
 l_back = [0,0,7.8,0,8,8.1,9.5,7.9,10.3,7.9];
 len_fleg = [0,0,19.9,0,19,19.3,21,19,17.5,19.3];
 len_bleg = [0,0,24.1,0,25,22.8,23,19.2,26.9,24.3];
 
% Swimming data

for g = 1:10 % this is the number of geckos we have (should stay constant)
	for t = 1:23 % this is the highest numbered trial we have for any gecko (for the purpose of reading in files)
		gecko_swim(g,t) = struct('F_slapstroke',{[]},'F_recoveryup',{[]},...
		'F_recoverydown',{[]},'B_slapstroke',{[]},'B_recoveryup',{[]},...
		'B_recoverydown',{[]},'Head',{[]}); % make a giant structure -- most of the entries will be blank if 
	end
end
waterlevel = zeros(num_trials*num_geckos,1);
for which_gecko = 1:num_geckos
	for which_trial = 1:num_trials
		fleg = strcat('gecko', num2str(which_gecko), '_s', num2str(which_trial), '_frontleg.txt');
		bleg = strcat('gecko', num2str(which_gecko), '_s', num2str(which_trial), '_backleg.txt');
		head = strcat('gecko', num2str(which_gecko), '_s', num2str(which_trial), '_head.txt');
		try	
			A = dlmread(fleg,',',11,0);
		catch err
			disp([fleg, ' does not exist']);
			continue
		end
		try	
			B = dlmread(bleg,',',11,0);
		catch err
			disp([bleg, ' does not exist']);
			continue
		end
		try	
			C = dlmread(head,',',11,0);
		catch err
			disp([head, ' does not exist']);
			continue
		end
		if isempty(find(A(:,4)~=-1,4))~=1 || isempty(find(B(:,4)~=-1,4))~=1 || isempty(find(C(:,4)~=-1,4))~=1
			waterlevel((which_gecko-1)*num_trials+which_trial) = mean(A(A(:,4)~=-1,4));
			currWL = waterlevel((which_gecko-1)*num_trials+which_trial);
		else
			disp(['gecko',num2str(which_gecko),'_trial',num2str(which_trial), ' did not have water level marked']);
			continue
		end
		
		indices = cell(7,1);
		fields = fieldnames(gecko_swim(which_gecko,which_trial));
		for i = 1:length(fieldnames(gecko_swim(which_gecko,which_trial)))
			start_column = 4;
			if i < 4
				indices{i} = find(A(:,i*2+start_column)~=-1);
				for k = 1:length(indices{i})
					gecko_swim(which_gecko,which_trial).(char(fields(i)))(end+1,:) = [A(indices{i}(k),1), A(indices{i}(k),i*2+start_column-1), abs(A(indices{i}(k),i*2+start_column)-currWL)];
				end
			elseif i < 7 && i > 3
				indices{i} = find(B(:,(i-3)*2+start_column)~=-1);
				for k = 1:length(indices{i})
					gecko_swim(which_gecko,which_trial).(char(fields(i)))(end+1,:) = [B(indices{i}(k),1), B(indices{i}(k),(i-3)*2+start_column-1), abs(B(indices{i}(k),(i-3)*2+start_column)-currWL)];
				end
			else
				indices{i} = find(C(:,(i-6)*2+start_column)~=-1);
				for k = 1:length(indices{i})
					gecko_swim(which_gecko,which_trial).(char(fields(i)))(end+1,:) = [C(indices{i}(k),1), C(indices{i}(k),(i-6)*2+start_column-1), abs(C(indices{i}(k),(i-6)*2+start_column)-currWL)];
				end
			end
		end
	end
end

%% calculate metrics we want
 
fwd_vel = cell(num_geckos,num_trials);
head_height = cell(num_geckos,num_trials);
frontlimb_height = cell(num_geckos,num_trials);
backlimb_height = cell(num_geckos,num_trials);
stance = cell(num_geckos,num_trials);
swing = cell(num_geckos,num_trials);
slap_vel = cell(num_geckos,num_trials);
stroke_vel_vert = cell(num_geckos,num_trials);
stroke_vel_hor = cell(num_geckos,num_trials);
bslap_vel = cell(num_geckos,num_trials);
bstroke_vel_vert = cell(num_geckos,num_trials);
bstroke_vel_hor = cell(num_geckos,num_trials);
period = cell(num_geckos,num_trials);



avg_vel = zeros(num_geckos,num_trials);
avgheadheight = zeros(num_geckos,num_trials);
avgFLheight = zeros(num_geckos,num_trials);
avgBLheight = zeros(num_geckos,num_trials);
avg_stance = zeros(num_geckos,num_trials);
avg_swing = zeros(num_geckos,num_trials);
avg_slapvel = zeros(num_geckos,num_trials);
avg_strokevel_vert = zeros(num_geckos,num_trials);
avg_strokevel_hor = zeros(num_geckos,num_trials);
avg_bslapvel = zeros(num_geckos,num_trials);
avg_bstrokevel_vert = zeros(num_geckos,num_trials);
avg_bstrokevel_hor = zeros(num_geckos,num_trials);
avg_period = zeros(num_geckos,num_trials);
fslap_imp = zeros(num_geckos,num_trials);
bslap_imp = zeros(num_geckos,num_trials);
fstr_imp_vert = zeros(num_geckos,num_trials);
fstr_imp_hor = zeros(num_geckos,num_trials);
bstr_imp_vert = zeros(num_geckos,num_trials);
bstr_imp_hor = zeros(num_geckos,num_trials);
min_imp = zeros(num_geckos,num_trials);
est_imp = zeros(num_geckos,num_trials);
percent = zeros(num_geckos,num_trials);

for which_gecko = 1:num_geckos
 	for which_trial = 1:num_trials
 		if isempty(gecko_swim(which_gecko,which_trial).Head) ~= 1
 			disp_x = diff(gecko_swim(which_gecko,which_trial).Head(:,2));
 			time = diff(gecko_swim(which_gecko,which_trial).Head(:,1))./500;
			
 			fwd_vel{which_gecko,which_trial} = abs(disp_x)./time;
			avg_vel(which_gecko,which_trial) = mean(fwd_vel{which_gecko,which_trial});
			
 			head_height{which_gecko,which_trial} = gecko_swim(which_gecko,which_trial).Head(:,3);
			avgheadheight(which_gecko,which_trial) = mean(head_height{which_gecko,which_trial});
			
 			frontlimb_height{which_gecko,which_trial} = gecko_swim(which_gecko,which_trial).F_recoverydown(:,3);
			avgFLheight(which_gecko,which_trial) = mean(frontlimb_height{which_gecko,which_trial});
			
 			backlimb_height{which_gecko,which_trial} = gecko_swim(which_gecko,which_trial).B_recoverydown(:,3);
			avgBLheight(which_gecko,which_trial) = mean(backlimb_height{which_gecko,which_trial});
			
			if gecko_swim(which_gecko,which_trial).F_recoveryup(1,1) > gecko_swim(which_gecko,which_trial).F_slapstroke(1,1)
				for i = 1:min(length(gecko_swim(which_gecko,which_trial).F_recoveryup(:,1)),length(gecko_swim(which_gecko,which_trial).F_slapstroke(:,1)))
					stance{which_gecko,which_trial}(end+1,:) = abs(gecko_swim(which_gecko,which_trial).F_recoveryup(i,1) - gecko_swim(which_gecko,which_trial).F_slapstroke(i,1))/500;
				end
			else
				for i = 1:min(length(gecko_swim(which_gecko,which_trial).F_recoveryup(:,1))-1,length(gecko_swim(which_gecko,which_trial).F_slapstroke(:,1)))
					stance{which_gecko,which_trial}(end+1,:) = abs(gecko_swim(which_gecko,which_trial).F_recoveryup(i+1,1) - gecko_swim(which_gecko,which_trial).F_slapstroke(i,1))/500;
				end
			end
			avg_stance(which_gecko,which_trial) = mean(stance{which_gecko,which_trial});
				
				
			if gecko_swim(which_gecko,which_trial).F_slapstroke(1,1) > gecko_swim(which_gecko,which_trial).F_recoveryup(1,1)
				for i = 1:min(length(gecko_swim(which_gecko,which_trial).F_slapstroke(:,1)),length(gecko_swim(which_gecko,which_trial).F_recoveryup(:,1)))
					swing{which_gecko,which_trial}(end+1,:) = abs(gecko_swim(which_gecko,which_trial).F_recoverydown(i,1) - gecko_swim(which_gecko,which_trial).F_slapstroke(i,1))/500;
				end
			else
				for i = 1:min(length(gecko_swim(which_gecko,which_trial).F_slapstroke(:,1))-1,length(gecko_swim(which_gecko,which_trial).F_recoveryup(:,1)))
					swing{which_gecko,which_trial}(end+1,:) = abs(gecko_swim(which_gecko,which_trial).F_slapstroke(i+1,1) - gecko_swim(which_gecko,which_trial).F_recoveryup(i,1))/500;
				end
			end
			avg_swing(which_gecko,which_trial) = mean(swing{which_gecko,which_trial});
			
			if gecko_swim(which_gecko,which_trial).F_slapstroke(1,1) > gecko_swim(which_gecko,which_trial).F_recoverydown(1,1)
				for i = 1:min(length(gecko_swim(which_gecko,which_trial).F_slapstroke(:,1)),length(gecko_swim(which_gecko,which_trial).F_recoverydown(:,1)))
					slap_vel{which_gecko,which_trial}(end+1,:) = abs(gecko_swim(which_gecko,which_trial).F_slapstroke(i,3) - gecko_swim(which_gecko,which_trial).F_recoverydown(i,3))/(abs(gecko_swim(which_gecko,which_trial).F_slapstroke(i,1) - gecko_swim(which_gecko,which_trial).F_recoverydown(i,1))/500);
				end
			else
				for i = 1:min(length(gecko_swim(which_gecko,which_trial).F_slapstroke(:,1))-1,length(gecko_swim(which_gecko,which_trial).F_recoverydown(:,1)))
					slap_vel{which_gecko,which_trial}(end+1,:) = abs(gecko_swim(which_gecko,which_trial).F_slapstroke(i+1,3) - gecko_swim(which_gecko,which_trial).F_recoverydown(i,3))/(abs(gecko_swim(which_gecko,which_trial).F_slapstroke(i+1,1) - gecko_swim(which_gecko,which_trial).F_recoverydown(i,1))/500);
				end
			end
			avg_slapvel(which_gecko,which_trial) = mean(slap_vel{which_gecko,which_trial});
	
			
			if gecko_swim(which_gecko,which_trial).F_recoveryup(1,1) > gecko_swim(which_gecko,which_trial).F_slapstroke(1,1)
				for i = 1:min(length(gecko_swim(which_gecko,which_trial).F_recoveryup(:,1)),length(gecko_swim(which_gecko,which_trial).F_slapstroke(:,1)))
					stroke_vel_vert{which_gecko,which_trial}(end+1,:) = abs(gecko_swim(which_gecko,which_trial).F_recoveryup(i,3) - gecko_swim(which_gecko,which_trial).F_slapstroke(i,3))/(abs(gecko_swim(which_gecko,which_trial).F_recoveryup(i,1) - gecko_swim(which_gecko,which_trial).F_slapstroke(i,1))/500);
                    stroke_vel_hor{which_gecko,which_trial}(end+1,:) = abs(gecko_swim(which_gecko,which_trial).F_recoveryup(i,2) - gecko_swim(which_gecko,which_trial).F_slapstroke(i,2))/(abs(gecko_swim(which_gecko,which_trial).F_recoveryup(i,1) - gecko_swim(which_gecko,which_trial).F_slapstroke(i,1))/500);

				end
			else
				for i = 1:min(length(gecko_swim(which_gecko,which_trial).F_recoveryup(:,1))-1,length(gecko_swim(which_gecko,which_trial).F_slapstroke(:,1)))
					stroke_vel_vert{which_gecko,which_trial}(end+1,:) = abs(gecko_swim(which_gecko,which_trial).F_recoveryup(i+1,3) - gecko_swim(which_gecko,which_trial).F_slapstroke(i,3))/(abs(gecko_swim(which_gecko,which_trial).F_recoveryup(i+1,1) - gecko_swim(which_gecko,which_trial).F_slapstroke(i,1))/500);
                    stroke_vel_hor{which_gecko,which_trial}(end+1,:) = abs(gecko_swim(which_gecko,which_trial).F_recoveryup(i+1,2) - gecko_swim(which_gecko,which_trial).F_slapstroke(i,2))/(abs(gecko_swim(which_gecko,which_trial).F_recoveryup(i+1,1) - gecko_swim(which_gecko,which_trial).F_slapstroke(i,1))/500);
				end
			end
			avg_strokevel_vert(which_gecko,which_trial) = mean(stroke_vel_vert{which_gecko,which_trial});
            avg_strokevel_hor(which_gecko,which_trial) = mean(stroke_vel_hor{which_gecko,which_trial});

			
			if gecko_swim(which_gecko,which_trial).B_recoveryup(1,1) > gecko_swim(which_gecko,which_trial).B_slapstroke(1,1)
				for i = 1:min(length(gecko_swim(which_gecko,which_trial).B_recoveryup(:,1)),length(gecko_swim(which_gecko,which_trial).B_slapstroke(:,1)))
					bstroke_vel_vert{which_gecko,which_trial}(end+1,:) = abs(gecko_swim(which_gecko,which_trial).B_recoveryup(i,3) - gecko_swim(which_gecko,which_trial).B_slapstroke(i,3))/(abs(gecko_swim(which_gecko,which_trial).B_recoveryup(i,1) - gecko_swim(which_gecko,which_trial).B_slapstroke(i,1))/500);
                    bstroke_vel_hor{which_gecko,which_trial}(end+1,:) = abs(gecko_swim(which_gecko,which_trial).B_recoveryup(i,2) - gecko_swim(which_gecko,which_trial).B_slapstroke(i,2))/(abs(gecko_swim(which_gecko,which_trial).B_recoveryup(i,1) - gecko_swim(which_gecko,which_trial).B_slapstroke(i,1))/500);

				end
			else
				for i = 1:min(length(gecko_swim(which_gecko,which_trial).B_recoveryup(:,1))-1,length(gecko_swim(which_gecko,which_trial).B_slapstroke(:,1)))
	               bstroke_vel_vert{which_gecko,which_trial}(end+1,:) = abs(gecko_swim(which_gecko,which_trial).B_recoveryup(i+1,3) - gecko_swim(which_gecko,which_trial).B_slapstroke(i,3))/(abs(gecko_swim(which_gecko,which_trial).B_recoveryup(i+1,1) - gecko_swim(which_gecko,which_trial).B_slapstroke(i,1))/500);
                   bstroke_vel_hor{which_gecko,which_trial}(end+1,:) = abs(gecko_swim(which_gecko,which_trial).B_recoveryup(i+1,2) - gecko_swim(which_gecko,which_trial).B_slapstroke(i,2))/(abs(gecko_swim(which_gecko,which_trial).B_recoveryup(i+1,1) - gecko_swim(which_gecko,which_trial).B_slapstroke(i,1))/500);
				end
			end
			avg_bstrokevel_vert(which_gecko,which_trial) = mean(bstroke_vel_vert{which_gecko,which_trial});
            avg_bstrokevel_hor(which_gecko,which_trial) = mean(bstroke_vel_hor{which_gecko,which_trial});
			
			if gecko_swim(which_gecko,which_trial).B_slapstroke(1,1) > gecko_swim(which_gecko,which_trial).B_recoverydown(1,1)
				for i = 1:min(length(gecko_swim(which_gecko,which_trial).B_slapstroke(:,1)),length(gecko_swim(which_gecko,which_trial).B_recoverydown(:,1)))
					bslap_vel{which_gecko,which_trial}(end+1,:) = abs(gecko_swim(which_gecko,which_trial).B_slapstroke(i,3) - gecko_swim(which_gecko,which_trial).B_recoverydown(i,3))/(abs(gecko_swim(which_gecko,which_trial).B_slapstroke(i,1) - gecko_swim(which_gecko,which_trial).B_recoverydown(i,1))/500);
				end
			else
				for i = 1:min(length(gecko_swim(which_gecko,which_trial).B_slapstroke(:,1))-1,length(gecko_swim(which_gecko,which_trial).B_recoveryup(:,1)))
					bslap_vel{which_gecko,which_trial}(end+1,:) = abs(gecko_swim(which_gecko,which_trial).B_slapstroke(i+1,3) - gecko_swim(which_gecko,which_trial).B_recoverydown(i,3))/(abs(gecko_swim(which_gecko,which_trial).B_slapstroke(i+1,1) - gecko_swim(which_gecko,which_trial).B_recoverydown(i,1))/500);
				end
			end
			avg_bslapvel(which_gecko,which_trial) = mean(bslap_vel{which_gecko,which_trial});
			
			
			period{which_gecko,which_trial} = diff(gecko_swim(which_gecko,which_trial).F_slapstroke(:,1))/500;
			avg_period(which_gecko,which_trial) = mean(period{which_gecko,which_trial});
			
			% impulse calculations
			C_d = 0.703; % drag coefficient -- look into changing later!
			farea = pi*(l_front./2000).^2; % area of front foot (as a disc) in m^2
			barea = pi*(l_back./2000).^2; % area of back foot (as a disc) in m^2
			
			fslap_imp(which_gecko,which_trial) = 4/3*999.97*(l_front(which_gecko)/2000)^3*mean(slap_vel{which_gecko,which_trial})/1000; % 4/3*rho*r^3*u_peak
			fstr_imp_vert(which_gecko,which_trial) = (1/2*C_d*farea(which_gecko)*999.97*(mean(stroke_vel_vert{which_gecko,which_trial})/1000)^2 + 1/2*C_d*farea(which_gecko)*999.97*9.81*len_fleg(which_gecko)/1000)*(len_fleg(which_gecko)/(mean(stroke_vel_vert{which_gecko,which_trial}))); % (0.5*C_d*S*rho*u_rms^2 + 0.5*C_d*S*rho*g*L)*(L/u_rms)
			fstr_imp_hor(which_gecko,which_trial) = (1/2*C_d*farea(which_gecko)*999.97*(mean(stroke_vel_hor{which_gecko,which_trial})/1000)^2 + 1/2*C_d*farea(which_gecko)*999.97*9.81*len_fleg(which_gecko)/1000)*(len_fleg(which_gecko)/(mean(stroke_vel_hor{which_gecko,which_trial}))); % (0.5*C_d*S*rho*u_rms^2 + 0.5*C_d*S*rho*g*L)*(L/u_rms)
            bslap_imp(which_gecko,which_trial) = 4/3*999.97*(l_back(which_gecko)/2000)^3*mean(bslap_vel{which_gecko,which_trial})/1000; % 4/3*rho*r^3*u_peak;
			bstr_imp_vert(which_gecko,which_trial) = (1/2*C_d*barea(which_gecko)*999.97*(mean(bstroke_vel_vert{which_gecko,which_trial})/1000)^2 + 1/2*C_d*barea(which_gecko)*999.97*9.81*len_bleg(which_gecko)/1000)*(len_bleg(which_gecko)/(mean(bstroke_vel_vert{which_gecko,which_trial}))); % (0.5*C_d*S*rho*u_rms^2 + 0.5*C_d*S*rho*g*L)*(L/u_rms)
			bstr_imp_hor(which_gecko,which_trial) = (1/2*C_d*barea(which_gecko)*999.97*(mean(bstroke_vel_hor{which_gecko,which_trial})/1000)^2 + 1/2*C_d*barea(which_gecko)*999.97*9.81*len_bleg(which_gecko)/1000)*(len_bleg(which_gecko)/(mean(bstroke_vel_hor{which_gecko,which_trial}))); % (0.5*C_d*S*rho*u_rms^2 + 0.5*C_d*S*rho*g*L)*(L/u_rms)
            min_imp(which_gecko,which_trial) = weights(which_gecko)/1000*9.81*avg_period(which_gecko,which_trial);
          	est_imp(which_gecko,which_trial) = fslap_imp(which_gecko,which_trial) + fstr_imp_vert(which_gecko,which_trial) + bslap_imp(which_gecko,which_trial) + bstr_imp_vert(which_gecko,which_trial);
			percent(which_gecko,which_trial) = est_imp(which_gecko,which_trial)/min_imp(which_gecko,which_trial);
			
 		end
 	end
end
%% output matrices in data files

dlmwrite('avg_vel.dat',avg_vel,'\t')
dlmwrite('avg_headheight.dat',avgheadheight,'\t')
dlmwrite('avg_FLheight.dat',avgFLheight,'\t')
dlmwrite('avg_BLheight.dat',avgBLheight,'\t')
dlmwrite('avg_stance.dat',avg_stance,'\t')
dlmwrite('avg_swing.dat',avg_swing,'\t')
dlmwrite('avg_fwdslapvel.dat',avg_slapvel,'\t')
dlmwrite('avg_fwdstrokevelvert.dat',avg_strokevel_vert,'\t')
dlmwrite('avg_fwdstrokevelhor.dat',avg_strokevel_hor,'\t')
dlmwrite('avg_backslapvel.dat',avg_bslapvel,'\t')
dlmwrite('avg_backstrokevelvert.dat',avg_bstrokevel_vert,'\t')
dlmwrite('avg_backstrokevelhor.dat',avg_bstrokevel_hor,'\t')
dlmwrite('avg_period.dat',avg_period,'\t')

%% put data into one big file
total_trials = length(find(avg_vel~=0));
compiled_data = zeros(total_trials,16);

i = 1; % index to count trials, so lazy -- redo in a nicer way, later.
for which_gecko = 1:num_geckos
	for which_trial = 1:num_trials
		if avg_vel(which_gecko,which_trial) ~= 0 
			compiled_data(i,1) = which_gecko;
			compiled_data(i,2) = which_trial;
			compiled_data(i,3) = avgFLheight(which_gecko,which_trial);
			compiled_data(i,4) = avgBLheight(which_gecko,which_trial);
			compiled_data(i,5) = avg_vel(which_gecko,which_trial);
			compiled_data(i,6) = avgheadheight(which_gecko,which_trial);
			compiled_data(i,7) = avg_slapvel(which_gecko,which_trial);
			compiled_data(i,8) = avg_bslapvel(which_gecko,which_trial);
			compiled_data(i,9) = avg_strokevel_vert(which_gecko,which_trial);
            compiled_data(i,10) = avg_strokevel_hor(which_gecko,which_trial);
			compiled_data(i,11) = avg_bstrokevel_vert(which_gecko,which_trial);
			compiled_data(i,12) = avg_bstrokevel_hor(which_gecko,which_trial);
			compiled_data(i,13) = avg_stance(which_gecko,which_trial);
			compiled_data(i,14) = avg_swing(which_gecko,which_trial);
			compiled_data(i,15) = avg_period(which_gecko,which_trial);
			compiled_data(i,16) = fslap_imp(which_gecko,which_trial);
			compiled_data(i,17) = fstr_imp_vert(which_gecko,which_trial);
            compiled_data(i,18) = fstr_imp_hor(which_gecko,which_trial);
			compiled_data(i,19) = bslap_imp(which_gecko,which_trial);
			compiled_data(i,20) = bstr_imp_vert(which_gecko,which_trial);
            compiled_data(i,21) = bstr_imp_hor(which_gecko,which_trial);
			compiled_data(i,22) = min_imp(which_gecko,which_trial);
			compiled_data(i,23) = est_imp(which_gecko,which_trial);
			compiled_data(i,24) = percent(which_gecko,which_trial);
			i = i+1;
		end
	end
end

dlmwrite('compiled_data.dat',compiled_data,'\t');
 
%% plotting
 
% figure;
% set(gca,'FontSize',18)
%  
% plot(avgheadheight,avg_vel,'*');
% P = polyfit(avgheadheight,avg_vel,1);
% hold on;
% plot(avgheadheight,P(1)*avgheadheight+P(2),'r','LineWidth',2);
% [rho,pval] = corr(avgheadheight',avg_vel');
% string1 = strcat('r = ', num2str(rho));
% string2 = strcat('p =  ', num2str(pval));
% text(20,950,string1,'FontSize',20,'FontName','Helvetica')
% text(20,900,string2,'FontSize',20,'FontName','Helvetica')
% xlabel('Head Height (mm)','FontSize',20,'FontName','Helvetica')
% ylabel('Forward Velocity (mm/s)','FontSize',20,'FontName','Helvetica')
% title('Head Height vs. Forward Velocity','FontSize',24,'FontName','Helvetica')
% 
% % Front limb height vs. forward velocity
%   
% figure;
% set(gca,'FontSize',18)
%  
% plot(avgFLheight,avg_vel,'*');
% P = polyfit(avgFLheight,avg_vel,1);
% hold on;
% plot(avgFLheight,P(1)*avgFLheight+P(2),'r','LineWidth',2);
% [rho,pval] = corr(avgFLheight',avg_vel');
% string1 = strcat('r = ', num2str(rho));
% string2 = strcat('p =  ', num2str(pval));
% text(15,950,string1,'FontSize',20,'FontName','Helvetica')
% text(15,900,string2,'FontSize',20,'FontName','Helvetica')
% xlabel('Front Limb Height (mm)','FontSize',20,'FontName','Helvetica')
% ylabel('Forward Velocity (mm/s)','FontSize',20,'FontName','Helvetica')
% title('Front Limb vs. Forward Velocity','FontSize',24,'FontName','Helvetica')
% 
% % Back limb height vs. front velocity
%   
% figure;
% set(gca,'FontSize',18)
% 
% plot(avgBLheight,avg_vel,'*');
% P = polyfit(avgBLheight,avg_vel,1);
% hold on;
% plot(avgBLheight,P(1)*avgBLheight+P(2),'r','LineWidth',2);
% [rho,pval] = corr(avgBLheight',avg_vel');
% string1 = strcat('r = ', num2str(rho));
% string2 = strcat('p =  ', num2str(pval));
% text(30,950,string1,'FontSize',20,'FontName','Helvetica')
% text(30,900,string2,'FontSize',20,'FontName','Helvetica')
% xlabel('Back Limb Height (mm)','FontSize',20,'FontName','Helvetica')
% ylabel('Forward Velocity (mm/s)','FontSize',20,'FontName','Helvetica')
% title('Back Limb Height vs. Forward Velocity','FontSize',24,'FontName','Helvetica')
% 
% %% fwd vel vs. duty factor in swimming
%  
% marks = {'FrontLeftEnter','FrontLeftExit','FrontRightEnter','FrontRightExit','BackLeftEnter',...
%  	'BackLeftExit','BackRightEnter','BackRightExit'};
%  
% df = zeros(length(marks)/4,27);
%  
% for which_gecko = 1:10
%  	for which_trial = 1:10
% 		for mark = 4:4:length(marks)
%  			enter = abs(gecko_DF_swim(which_gecko,which_trial).(char(marks(mark-1))));
% 			exit = abs(gecko_DF_swim(which_gecko,which_trial).(char(marks(mark))));
% 			if isempty(gecko_swim(which_gecko,which_trial).side.Head) == 1
% 				df(:,(which_gecko-1)*3+which_trial) = 100;
% 				continue
% 			end
% 			if isempty(exit) == 1 || isempty(enter) == 1
% 				continue
% 			end
% 			if enter(1,1) < exit(1,1)
% 				first_exit = 2;
% 			else
% 				first_exit = 1;
% 			end
% 			temp = [];
% 			for i = 1:min(length(exit(:,1))-first_exit+1,length(enter(:,1))-1)
% 				temp(i) = abs((enter(i,1) - exit(first_exit+i-1,1))/(enter(i+1,1)-enter(i,1)));
% 			end
% 			df(ceil(mark/4),(which_gecko-1)*3+which_trial) = mean(temp);
% 		end
% 	end
% end
% 
% dff_full = df(1,df(1,:)~=100);
% dfb_full = df(2,df(2,:)~=100);
% 
% dff = dff_full(dff_full > 0);
% dfb = dfb_full(dff_full > 0);
% 
% figure;
% set(gca,'FontSize',18)
% 
% plot(dff,avg_vel,'*')
% P = polyfit(dff,avg_vel,1);
% hold on;
% plot(dff,P(1)*dff+P(2),'r','LineWidth',2);
% [rho,pval] = corr(dff',avg_vel');
% string1 = strcat('r = ', num2str(rho));
% string2 = strcat('p =  ', num2str(pval));
% text(0.6,800,string1,'FontSize',20,'FontName','Helvetica')
% text(0.6,750,string2,'FontSize',20,'FontName','Helvetica')
% xlabel('Duty Factor', 'FontSize',20,'FontName','Helvetica')
% ylabel('Forward Velocity (mm/s)', 'FontSize',20,'FontName','Helvetica')
% title('Duty Factor of Front Limb vs. Forward Velocity','FontSize',24)
%  
% figure;
% set(gca,'FontSize',18)
% 
% plot(dfb,avg_vel,'*')
% P = polyfit(dfb,avg_vel,1);
% hold on;
% plot(dfb,P(1)*dfb+P(2),'r','LineWidth',2);
% [rho,pval] = corr(dfb',avg_vel_swim');
% string1 = strcat('r = ', num2str(rho));
% string2 = strcat('p =  ', num2str(pval));
% text(0.6,800,string1,'FontSize',20,'FontName','Helvetica')
% text(0.6,750,string2,'FontSize',20,'FontName','Helvetica')
% xlabel('Duty Factor', 'FontSize',20,'FontName','Helvetica')
% ylabel('Forward Velocity (mm/s)', 'FontSize',20,'FontName','Helvetica')
% title('Duty Factor of Back Limb vs. Forward Velocity','FontSize',24)
