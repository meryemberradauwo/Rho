clc;clear;close all;
dbstop if error
% Reading Data
file_name = 'exampledataset.xlsx';
data = importdata(file_name);
data = transpose(data.data);

%thresholds to be modified for each file (n2 is a more complexe variable defined later)
n=2; %2  relates to n2
n1=-1e-4; %-1e-4 relates to amplitude of negative voltage values
n3 = 10; %20  controls outliers, max difference between T and V values CHANGE THIS **********
n4=1.3e-4; %1.3e-4 differentiates TB, TA from difference between T and V values
n5=3; %3 minimum number of elements to consider
n6=10; %number of elements for T, dont change

%constants
I=0.2; % constant current [A]
Lnum=2.44*10^(-8); % Lorenz number [W*Ohm*K^-2]
L=0.001620; % sample length  [m] (60t) 171
D=0.000508; % sample diameter [m] (20t) 508
errL=0.0000127; % error from microscope [m]
errD=errL; 
Lp=0.00001017; % disc length [m] (4t)
Dp=0.00127; % disc diameter [m] (50t)

% plot of raw data
figure(1);
plot(data,'.b'); 
title('Step1 : Temporary result');
xlabel('Index')
ylabel('Data')
grid on;grid minor;

%This identifies all negative values
a = data;
b = [data(2:end) 0];
c = [data(3:end) 0 0];
d = [data(4:end) 0 0 0];
ab = a.*b; 
ac = a.*c;
ad = a.*d;
negative_index = ab<0 & ac<0 & ad<0; % finds if ab, ac and ad all contain negative numbers at same locations
negative_index_number = (find(negative_index)); %finds the values that are negatives
aa = a(negative_index_number);

%selecting sections based on finding the sections with negative numbers
k = 1;
neg_diff = []; %create an empty matrix
interval_neg2=[]; %create an empty matrix

%loop to identify when the first voltage is positive or negative
for i=1:(length(negative_index_number)-1)
        interval_neg = (negative_index_number(i)+1:negative_index_number(i+1));
        if sum(data(interval_neg)<0)>n5 && mode(data(interval_neg))<n1 %mode returns the most frequent value
            temp_data_neg = data(interval_neg); %defines a vector of negative values
            temp_data_neg = temp_data_neg(temp_data_neg<0);
            interval_neg = interval_neg(temp_data_neg<0); % inditifies all sections of negative values
            
            %updates figure 1 by highlighting negative values in red
            figure(1);
            hold on;
            plot(interval_neg,temp_data_neg,'.r');
            
            temp_data_before = data(interval_neg(1)-n6:interval_neg(1)-1); %selects n6 values before first negative values
            temp_data_after = data(interval_neg(end)+1:interval_neg(end)+n6); %selects n6 values after first negative values
            
            if k==44
                1;
            end
            
     % in the case where the sequence is TB VP VN TA
            
            %If the difference between Tb-negative values with Ta-negative 
            %values is less than n4, then interval of negative values is
            %defined from first to last element of negative values
            %whereas if Tb-negative values is smaller than Ta-negative 
            %values, then the positive values are first part of negative values
            
            n2 = max(ceil(n*length(interval_neg)),n6); %ceil rounds each entry of bracket, max finds the highest number
            if abs(abs(abs(median(temp_data_before))-abs(median(temp_data_neg)))- abs(abs(median(temp_data_after))-abs(median(temp_data_neg))))<n4
                interval_neg2 = [interval_neg2 interval_neg(1) interval_neg(end)];
            elseif (abs(abs(median(temp_data_before))-abs(median(temp_data_neg))) < abs(abs(median(temp_data_after))-abs(median(temp_data_neg))))
                interval_pos = (interval_neg(1)-n2:interval_neg(1)-1); %first part of negative values
                temp_data_pos = data(interval_pos); %defines positive values as y axis values
                V_ind = (min([interval_pos , interval_neg]):max([interval_pos , interval_neg])); %creates one vector of all voltage values
                VP_ind = V_ind(data(V_ind)>=0); %defines positive voltage
                isempty(VP_ind); %find empty cells, should be 0  
                VP = data(VP_ind); %VP is y values of the VP indices 
                VN_ind = V_ind(data(V_ind)<0); %defines negative voltage
                VN = data(VN_ind); %VN is y values of the VN indices
                dif_VN = diff(VN); %calculates difference between adjacent elements in VN
                T_a_ind = (VN_ind(end)+1:VN_ind(end)+n6); %defines indexes for Ta as n6 values after VN
                T_a = data(T_a_ind); %T_a is y values of Ta indices
                
                %repeat the previous for each data point
                [~,ind1] = max(abs(diff(VP(1:end-n5))));
                if ~isempty(ind1)
                    ind = VP_ind(1)+ind1;
                    VP_ind = (ind:VP_ind(end));
                    VP = data(VP_ind);
                    T_b_ind = ind-n6:ind-1;
                    T_b = data(T_b_ind);
                    V_ind = (VP_ind(1):VN_ind(end)); %creates a vector of all VP and VN data points
                    final.VP{k} = [VP;VP_ind]; %creates vector of all VP data points
                    final.VN{k} = [VN;VN_ind]; %creates vector of all VN data points
                    final.Tb{k} = [T_b;T_b_ind]; %creates vector of all Tb data points
                    final.Ta{k} = [T_a;T_a_ind]; %creates vector of all Ta data points
                    k = k+1;
                    
                    %updates figure 1 to identify all data points with this sentence
                    plot(V_ind,data(V_ind),'*k',T_b_ind,T_b,'.r',T_a_ind,T_a,'.g');
                end
                
     % in the case where the sequence is TB VN VP TA
            
            %define positive values as second part of negative values
            else 
                interval_pos = (interval_neg(end)+1:interval_neg(end)+n2); %second part of negative values
                temp_data_pos = data(interval_pos); %defines positive values as y axis values
                V_ind = (min([interval_pos , interval_neg]):max([interval_pos , interval_neg]));%creates one vector of all voltage values
                VP_ind = V_ind(data(V_ind)>=0);%defines positive voltage
                isempty(VP_ind); %find empty cells, should be 0
                VP = data(VP_ind); %VP is y values of the VP indices 
                VN_ind = V_ind(data(V_ind)<0); %defines negative voltage
                VN = data(VN_ind); %VN is y values of the VN indices
                dif_VN = diff(VN); %calculates difference between adjacent elements in VN
                if ~isempty(VN_ind)
                    T_b_ind = (VN_ind(1)-n6:VN_ind(1)-1); %defines Tb as n6 values before VN
                    T_b = data(T_b_ind); %defines y values of Tb
                    
                     %repeat the previous for each data point
                    [~,ind1] = max(abs(diff(VP(n5:end))));
                    if ~isempty(ind1)
                        ind = VP_ind(1)+ind1+(n5-1);
                        VP_ind = (VP_ind(1):ind-1);
                        VP = data(VP_ind);
                        T_a_ind = ind:ind+(n6-1);
                        T_a = data(T_a_ind);
                        V_ind = (VN_ind(1):VP_ind(end));
                        final.VP{k} = [VP;VP_ind];
                        final.VN{k} = [VN;VN_ind];
                        final.Tb{k} = [T_b;T_b_ind];
                        final.Ta{k} = [T_a;T_a_ind];
                        k = k+1;
                        
                        %updates figure 1 to identify all data points with this sentence
                        plot(V_ind,data(V_ind),'*k',T_b_ind,T_b,'.r',T_a_ind,T_a,'.g');
                    end
                end
            end            
        end
end

%loop to identify Ta and Tb depending on the previous scenarios
for jj=1:length(interval_neg2)-1
   interval_neg = interval_neg2(jj):interval_neg2(jj+1);
   n2 = max(ceil(n*length(interval_neg)),n6);
   if length(interval_neg)>1 && mode(data(interval_neg))<n1
        temp_data_neg = data(interval_neg); %defines a vector of negative values
        temp_data_neg = temp_data_neg(temp_data_neg<0);
        interval_neg = interval_neg(temp_data_neg<0); %identifies all sections of negative values
        
        %updates figure 1 by highlighting negative values in red
        figure(1);
        hold on;
        plot(interval_neg,temp_data_neg,'.r')
        
        temp_data_before = data(interval_neg(1)-n2:interval_neg(1)-1); %selects first part of negative values
        temp_data_after = data(interval_neg(end)+1:interval_neg(end)+n2); %selects second part of negative values
        
     % in the case where the sequence is TB VP VN TA
     
        %if Tb>Ta then positive interval is first part of negative values
        if max(abs(diff(temp_data_before)))>max(abs(diff(temp_data_after)))
            interval_pos = (interval_neg(1)-n2:interval_neg(1)-1); %first part of negative values
            temp_data_pos = data(interval_pos); %defines y axis values of positive values
            V_ind = (min([interval_pos , interval_neg]):max([interval_pos , interval_neg])); %creates one vector of all negative values
            VP_ind = V_ind(data(V_ind)>=0); %defines positive voltage
            isempty(VP_ind); %find empty cells, should be 0
            VP = data(VP_ind); %VP is y values of the VP indices
            VN_ind = V_ind(data(V_ind)<0); %defines negative voltages
            VN = data(VN_ind); %VN is y values of the VN indices
            dif_VN = diff(VN); %calculate differences between adjacent elements in VN
            T_a_ind = (VN_ind(end)+1:VN_ind(end)+n6); %defines Ta as n6 values after VN
            T_a = data(T_a_ind); %defines y values of Ta
            
            %repeat the previous for each data point
            [~,ind1] = max(abs(diff(VP(1:end-n5))));
            if ~isempty(ind1)
                ind = VP_ind(1)+ind1;
                VP_ind = (ind:VP_ind(end));
                VP = data(VP_ind);
                T_b_ind = ind-n6:ind-1;
                T_b = data(T_b_ind);
                V_ind = (VP_ind(1):VN_ind(end));
                final.VP{k} = [VP;VP_ind];
                final.VN{k} = [VN;VN_ind];
                final.Tb{k} = [T_b;T_b_ind];
                final.Ta{k} = [T_a;T_a_ind];
                k = k+1;
                
                %updates figure 1 to identify all data points with this sentence
                plot(V_ind,data(V_ind),'*k',T_b_ind,T_b,'.r',T_a_ind,T_a,'.g');
            end
            
     % in the case where the sequence is TB VN VP TA
        else
            interval_pos = (interval_neg(end)+1:interval_neg(end)+n2); %second part of negative values
            temp_data_pos = data(interval_pos); %defines positive values as y axis values
            V_ind = (min([interval_pos , interval_neg]):max([interval_pos , interval_neg])); %creates one vector of all negative values
            VP_ind = V_ind(data(V_ind)>=0); %defines positive voltage
            isempty(VP_ind); %find empty cells, should be 0
            VP = data(VP_ind); %VP is y values of the VP indices
            VN_ind = V_ind(data(V_ind)<0); %defines negative voltage
            VN = data(VN_ind); %VN is y values of the VN indices
            dif_VN = diff(VN); %calculates differences between adjacent elements in VN
            if ~isempty(VN_ind)
                T_b_ind = (VN_ind(1)-n6:VN_ind(1)-1); %defines Tb as n6 values before VN
                T_b = data(T_b_ind); %defines y values of Tb
                
                %repeat the previous for each data point
                [~,ind1] = max(abs(diff(VP(n5:end)))); 
                if ~isempty(ind1)
                    ind = VP_ind(1)+ind1+(n5-1);
                    VP_ind = (VP_ind(1):ind-1);
                    VP = data(VP_ind);
                    T_a_ind = ind:ind+(n6-1);
                    T_a = data(T_a_ind);
                    V_ind = (VN_ind(1):VP_ind(end));
                    final.VP{k} = [VP;VP_ind];
                    final.VN{k} = [VN;VN_ind];
                    final.Tb{k} = [T_b;T_b_ind];
                    final.Ta{k} = [T_a;T_a_ind];
                    k = k+1;
                    
                    %updates figure 1 to identify all data points with this sentence
                    plot(V_ind,data(V_ind),'*k',T_b_ind,T_b,'.r',T_a_ind,T_a,'.g');
                end
            end
        end
   end
end

%plot of final selections, without outliers 
figure(2);
plot(data,'.b'); hold on;
title('Step2 : Final Result')
xlabel('Index')
ylabel('Data')
grid on;grid minor;

VP_dif_max = []; %create an empty matrix
VN_dif_max = []; %create an empty matrix
Ta_dif_max = []; %create an empty matrix
Tb_dif_max = []; %create an empty matrix

%loop for creating final variables without outliers
for kk = 1:length(final.VP)
     if kk==1
         Ta_ind = final.Ta{kk}(2,:); %defines TA
         VP_b_ind = final.VP{kk+1}(2,:); %scenario TB VP VN TA
         VN_b_ind = final.VN{kk+1}(2,:); %scenario TB VN VP TA
         for i=1:length(Ta_ind)
             item = Ta_ind(i);
             %if VP or VN indices overlap with TA, then indice is -1 
             if any(VP_b_ind==item) || any(VN_b_ind==item) 
                 Ta_ind(i) = -1;
             end
         end
         %all indices not equal to -1 are labelled as TA
         Ta = [final.Ta{kk}(1,Ta_ind~=-1);Ta_ind(Ta_ind~=-1)];
         final.Ta{kk} = Ta;
         
     elseif kk==length(final.VP) %for the rest of VP elements
         Tb_ind = final.Tb{kk}(2,:); %defines TB
         VP_a_ind = final.VP{kk-1}(2,:); %scenario TB VP VN TA
         VN_a_ind = final.VN{kk-1}(2,:); %scenario TB VN VP TA
         for i=1:length(Ta_ind)
             item = Tb_ind(i);
             %if VP or VN indices overlap with TB, then indice is -1
             if any(VP_a_ind==item) || any(VN_a_ind==item)
                 Tb_ind(i) = -1;
             end
         end
         %all indices not equal to -1 are labelled as TB
         Tb = [final.Tb{kk}(1,Tb_ind~=-1);Tb_ind(Tb_ind~=-1)];
         final.Tb{kk} = Tb;
         
     else 
        Ta_ind = final.Ta{kk}(2,:); %opposite scenario as previous loop for TA
        VP_a_ind = final.VP{kk-1}(2,:);
        VN_a_ind = final.VN{kk-1}(2,:);
        
        Tb_ind = final.Tb{kk}(2,:); %opposite scenario as previous loop for TB
        VP_b_ind = final.VP{kk+1}(2,:);
        VN_b_ind = final.VN{kk+1}(2,:);
                       
        for i=1:length(Ta_ind) %same as previous loop but for opposite scenario
           item = Ta_ind(i);
           if any(VP_b_ind==item) || any(VP_a_ind==item) || any(VN_b_ind==item) || any(VN_a_ind==item)
               Ta_ind(i) = -1;
           end
                          
           item = Tb_ind(i); %same as previous loop but for opposite scenario
           if any(VP_b_ind==item) || any(VP_a_ind==item) || any(VN_b_ind==item) || any(VN_a_ind==item)
               Tb_ind(i) = -1;
           end
        end
        
        %rewrites all TA and TB by avoiding overlapping indices
        Ta = [final.Ta{kk}(1,Ta_ind~=-1);Ta_ind(Ta_ind~=-1)];
        final.Ta{kk} = Ta;
        Tb = [final.Tb{kk}(1,Tb_ind~=-1);Tb_ind(Tb_ind~=-1)];
        final.Tb{kk} = Tb;
     end
     
     % identify y-axis values as first row of matrix final
     VP = final.VP{kk}(1,:); 
     VN = final.VN{kk}(1,:); 
     Ta = final.Ta{kk}(1,:); 
     Tb = final.Tb{kk}(1,:); 
     
     % identify x-axis values as second row of matrix final
     VP_ind = final.VP{kk}(2,:);
     VN_ind = final.VN{kk}(2,:);
     Ta_ind = final.Ta{kk}(2,:);
     Tb_ind = final.Tb{kk}(2,:);

     % rewrite VP without including values that have a difference larger than n3
     VP_dif = sort(abs(diff(VP)));
     if length(VP_dif)>n5
         VP_dif = VP_dif(end-n5);
         VP_dif1 = [abs(diff(VP)) max(abs(diff(VP)))];
         VP_dif2 = [max(abs(diff(VP))) abs(diff(VP))];
         VP_del_ind = VP_dif1>n3*VP_dif & VP_dif2>n3*VP_dif; 
         VP = VP(~VP_del_ind); 
         VP_ind = VP_ind(~VP_del_ind);
     end
     
     % rewrite VN without including values that have a difference larger than n3
     VN_dif = sort(abs(diff(VN)));
     if length(VN_dif)>n5
         VN_dif = VN_dif(end-n5);
         VN_dif1 = [abs(diff(VN)) max(abs(diff(VN)))];
         VN_dif2 = [max(abs(diff(VN))) abs(diff(VN))];
         VN_del_ind = VN_dif1>n3*VN_dif & VN_dif2>n3*VN_dif;
         VN = VN(~VN_del_ind);
         VN_ind = VN_ind(~VN_del_ind);
     end
     
     % rewrite TA without including values that have a difference larger than n3
     Ta_dif = sort(abs(diff(Ta)));
     if length(Ta_dif)>n5
         Ta_dif = Ta_dif(end-n5);
         Ta_dif1 = [abs(diff(Ta)) max(abs(diff(Ta)))];
         Ta_dif2 = [max(abs(diff(Ta))) abs(diff(Ta))];
         Ta_del_ind = Ta_dif1>n3*Ta_dif & Ta_dif2>n3*Ta_dif;
         Ta = Ta(~Ta_del_ind);
         Ta_ind = Ta_ind(~Ta_del_ind);
     end
     
     % rewrite TB without including values that have a difference larger than n3
     Tb_dif = sort(abs(diff(Tb)));
     if length(Tb_dif)>n5
         Tb_dif = Tb_dif(end-n5);
         Tb_dif1 = [abs(diff(Tb)) max(abs(diff(Tb)))];
         Tb_dif2 = [max(abs(diff(Tb))) abs(diff(Tb))];
         Tb_del_ind = Tb_dif1>n3*Tb_dif & Tb_dif2>n3*Tb_dif;
         Tb = Tb(~Tb_del_ind);
         Tb_ind = Tb_ind(~Tb_del_ind);
     end
          
     % rewrite final matrix
     final.VP{kk} = [VP;VP_ind];
     final.VN{kk} = [VN;VN_ind];
     final.Ta{kk} = [Ta;Ta_ind];
     final.Tb{kk} = [Tb;Tb_ind];
    
     % calculating mean values for each data point
     final_avg.VP(kk) = mean(VP);
     final_avg.VN(kk) = mean(VN);
     final_avg.Ta(kk) = mean(Ta);
     final_avg.Tb(kk) = mean(Tb);
     
     % calculating standard deviation for each data point
     final_avg.errVP(kk) = std(VP);
     final_avg.errVN(kk) = std(VN);
     final_avg.errTa(kk) = std(Ta);
     final_avg.errTb(kk) = std(Tb);
     
     % update figure 2 using final variables
     plot(final.VP{kk}(2,:),final.VP{kk}(1,:),'*k',...
         final.VN{kk}(2,:),final.VN{kk}(1,:),'*r',...
         final.Ta{kk}(2,:),final.Ta{kk}(1,:),'.r',...
         final.Tb{kk}(2,:),final.Tb{kk}(1,:),'.g');
     
end

%creating the final variables that will be used in equations
Vpositive=final_avg.VP;
Vnegative=final_avg.VN;
Tbefore=final_avg.Tb;
Tafter=final_avg.Ta;

%replacing outliers by interpolation
Vpos=filloutliers(Vpositive,'linear','mean');
Vneg=filloutliers(Vnegative,'linear','mean');
Tbef=filloutliers(Tbefore,'linear','mean');
Taft=filloutliers(Tafter,'linear','mean');

% Convert emf to temperature
emf=1/2.*(Taft+Tbef); % vector of emf value for each data point
Td=emf.*1000+0.342;
Tc=74.124732.*Td-4.28082813.*Td.^2+0.52113892.*Td.^3 ...
    -0.0457487201.*Td.^4+0.00280578284.*Td.^5 ...
    -0.000113145137.*Td.^6+0.00000285489684.*Td.^7 ...
    -0.0000000407643828.*Td.^8+0.000000000251358071.*Td.^9; 
T=Tc+273.15; % gives T in kelvin

Vt=1/2.*abs(Vpos-Vneg); % total voltage for each data point 
pW=2*10^(-6).*T.^2+0.0086.*T-0.4985; % fit to Josh's W 5Gpa data
Va=pW*10^-8*Lp*I*4/(pi*Dp^2); % V contribution from 1 plug
Vfinal=Vt-2*Va; % V from sample 
p=10^8*(pi.*D.^2.*Vfinal)./(4.*L.*I); % resistivity [microOhm*cm]
k=10^8*Lnum.*T./p; % thermal conductivity [W*K^-1*m^-1]

%error bar calculations
errVP=final_avg.errVP;
errVN=final_avg.errVN;
errTa=final_avg.errTa;
errTb=final_avg.errTb;

% Error Propagation, assumes no error on Lorenz number or current
errT=1/2.*(errTa+errTb); % error on T
errV=1/2*sqrt(errVP.^2+errVN.^2); % error on V 
errp=p.*sqrt((2*errD./D).^2+(errV./Vfinal).^2+(errL./L).^2); % error p in cm 
errk=errp/p*k; % error on k

%figure of resistivity data
figure (3)
plot(T,p,'*') % x axis is temperature 
xlabel('Temperature (K)')
ylabel('Resistivity (microO*cm)')
