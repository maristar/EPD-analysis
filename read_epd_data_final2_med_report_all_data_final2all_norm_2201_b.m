%% EPD data analysis
% This skript reads csv files from old types EPD in Rikshospitalet 
% and plots them. Then adds the figures in a report. 
% This opens also data from txt and plots them. And adds them in the same
%report.  Maria Stavrinou, Desember 2020, Rikshospitalet. 
% Revised 15 Januar 2021. 

clear all 
close all;

% Excels &  data location 
global_path ='C:\Users\Maria\Documents\Min_NUK_på_Asus\EPD_presentasjon\All_datas\'; 
global_raw_path='C:\Users\Maria\Documents\Min_NUK_på_Asus\EPD_presentasjon\All_datas\raw_data\';
cd(global_path);
mkdir('Figures_all_2021');
mkdir('Report_all_2021');
global_figures_path=[global_path 'Figures_all_2021\'];
global_report_path=[global_path 'Report_all_2021\'];
global_B = datetime(2020,01,01);
string_xlabel=[];

% Start the report
% TODO: To check if report can be saved in another directory - now it is saved in
% the raw data directory (global_raw_path).
import mlreportgen.report.*;
import mlreportgen.dom.*;
rpt=Report('EPD_report_all_2021_testmars21', 'pdf');
open(rpt)

%% First page 
T=Text("Appendix: EPD rapport 01 Januar 2020 - 27 november 2020.");
T.Bold=true;
T.FontSize="28pt";
add(rpt, T);

%% 1. Sjekk de gamle EPDS som har *.csv extension
Folder=cd(global_raw_path);
global_list = dir('*.csv');
global_N_File=numel(global_list);

for o=1:global_N_File %Loop start 
    qs=numel(global_list)-o; % where do i need that`?
    %clc%
    cd(global_raw_path);
    excelfile_name=global_list(o).name;
    disp(['Working on: ' global_list(o).name])
    ansatte_navn=global_list(o).name(1:end-4);
    % Get some info on imported file
    % opts = detectImportOptions(excelfile_name);
    
    %% Create chapters for each employee
    ch=Chapter;
    ch.Title=ansatte_navn;
    ch.Numbered=true; 
    add(rpt, ch);
    
  %% Make a table from the csv file
  % like in excel - add 2nd +0.0013rd column for Hp10 og 4th+0.0001*5th column 
    T2 = readtable(excelfile_name);
    %readtable('AlexanderSherwaniSave As.csv');
    SizeTable= size(T2);
    rows=SizeTable(1);
    
    for i=1:rows
    tres2(:,i)=T2.Var2(i)+(T2.Var3(i))*0.001;
    end
    clear i     
    % Create a new column with 
    T2.Hp10 =tres2';
    clear tres2;
   
    % Do the same with Hp07
    for i=1:rows
        temp_data2(:,i)=T2.Var4(i)+(T2.Var5(i))*0.001;
    end
    T2.Hp07 =temp_data2';
    clear temp_data2 i;
    
    % Make a plotable time axis
    if o==5, T2.Var1{408} = '15.08.2020 00:00:00'; end % Susanne hadde en rar timepoint
   for j=1:rows
        temp_time(:,j)=datetime(T2.Var1(j), 'InputFormat', 'dd.MM.yyyy HH:mm:ss');
    end
    time_vec=temp_time';
    clear j 
    clear temp_time; 
    % Cut out dates before 2020. 
    
    %% Ny av 15.01.2021 Define number of days 
    v=datevec(time_vec);
   % take the dates
    var10=datetime(v(:,1:3));
    % take the times
    var11=duration(v(:,4:end));
    % write them in the final table T2
    T2.var8=var10; T2.var9=var11; T2.timevec=time_vec;
   % find the unique indices of the 
    [C, ic]=unique(var10);
    number_of_days=length(C);
    T2.total_days(1,1)=number_of_days;
    %% find total number of days above 01.01.2020
    idx = find(time_vec>global_B);
    v2020=datevec(time_vec(idx));
    Alldays_2020=datetime(v2020(:,1:3));
    [C2020, ic2020]=unique(Alldays_2020);
    number_of_2020days=length(C2020);
    T2.total_days2020(1,1)=number_of_2020days;
    % end
    fprintf([global_list(o).name ': total days:' num2str(number_of_days) '. Total days in 2020:' num2str(number_of_2020days)]);
    %% Dates above 01.01.2020
    % New Friday 22 jan 2021
    
    idx = find(time_vec>global_B & T2.Hp10>0 &T2.Hp07>0);
    Data_to_plot_Hp10(:,:) = T2.Hp10(idx); % Pull out the data based on the idx
    Data_to_plot_Hp07(:,:) = T2.Hp07(idx); 
    time_vec_to_plot(:,:) =time_vec(idx);
   
   %% Normalize extracting the initial value of 01.01.2020 (new 18012020)
   firstHp10_2020=T2.Hp10(idx(1));
   firstHp07_2020=T2.Hp07(idx(1));
   
   size_norm=length(idx);
   
   % make a normalized column
   for kk=1:size_norm;
       index_norm=idx(1)+kk-1;
        T2.Hp10norm(index_norm)=T2.Hp10(index_norm)-firstHp10_2020;
        T2.Hp07norm(index_norm)=T2.Hp07(index_norm)-firstHp07_2020;
   end
    clear kk index_norm
  % Fill the other cells from 0 to idx-1 with zeros
  for jj=1:(idx(1)-1)
      T2.Hp10norm(jj)=0;
      T2.Hp07norm(jj)=0;
  end
  clear jj
  
  % Create Data_to_plot_norm variables 
  Data_to_plot_Hp10_norm(:,:)=T2.Hp10norm(idx);
  Data_to_plot_Hp07_norm(:,:)=T2.Hp07norm(idx);
  
  
  % Siste data: 
  Siste_Hp07=Data_to_plot_Hp07_norm(end);
  Siste_Hp10=Data_to_plot_Hp10_norm(end);
  
  SumDoseMatlab.(ansatte_navn).SisteHp07=Siste_Hp07;
  SumDoseMatlab.(ansatte_navn).SisteHp10=Siste_Hp10;
  
  SumDoseMatlab.(ansatte_navn).totaldays=number_of_2020days;
  
  T2.DoseSumHp07_sistenorm(1,1)=Siste_Hp07;
  T2.DoseSumHp10_sistenorm(1,1)=Siste_Hp10;
  
  % Check that first value of normalized values of idx(1) is zero 
  
  % Now clear idx
  clear idx size_norm firstHp10_2020 first_Hp07_2020
  % End normalize
    
  %% FIGURE 
    fig=figure; 
    plot(time_vec_to_plot, Data_to_plot_Hp10_norm, 'b*'); hold on; plot(time_vec_to_plot,Data_to_plot_Hp07_norm, 'r.');
    clear time_vec_to_plot
    %% For the figures
    % Make figures better
     legend('Hp10','Hp07','Location','northwest'); 
     xlabel('time');
     ylabel('Hp10(mSv) og Hp07(mSv)');
      ylim([0 6]);
    % xlim([01-Jan-2020 05:00:00 27-Nov-2020 20:00:00])
     time1=datetime("01.01.2020 00:00:09",'InputFormat', "dd.MM.yyyy HH:mm:ss");
     time2=datetime("27.11.2020 14:00:00", 'InputFormat', "dd.MM.yyyy HH:mm:ss");
     xlim([time1 time2]); 
     % Figure annotation which does not work 
     str_text={[global_list(o).name ' Total days in 2020:' num2str(number_of_2020days)]}
     text(2,7,str_text);
     %% For anonymt plot
     %title(['Hp10 ansatt nr '  num2str(o)]);
     
     %% For ikke anonymt plot
     title([ansatte_navn]);
     
     %% Save figures
     savename=[excelfile_name(1:end-5) '_results'];
     cd(global_path);
     cd(global_figures_path)
     mkdir(savename)
     cd(savename)
     savefig(fig,savename)
     saveas(fig, savename, 'tiff')
     saveas(fig, savename, 'fig')
     
     %% Add figure to rapport 
     fig1=Figure(fig);
     fig1.Width="4.2in";
     fig1.Height="4in";
     %fig1.Scaling="custom";
     % Add fig to rapport
     add(rpt, fig1); 
     close(fig)
     
     %% Sum Dose figures
     fig100=figure(100); plot(Siste_Hp07, 'r*'); hold on; plot(Siste_Hp10, 'b*'); hold on;
     ylim([0 20]); 
     
     
%      %string_xlabel=[string_xlabel (ansatte_navn)]
%      xticks('auto');
%      if o==1
%          string_labels=[];
%      end
%      
%      string_labels=[(string_labels), (ansatte_navn)]
%      xticklabels({(string_labels)}

     save T2_gamleEPD T2
     clear SizeTable rows Data_to_plot_Hp07 Data_to_plot_Hp10 idx time_vec time_vec_to_plot T2 excelfile_name savename qs rows columns
     clear Data_to_plot_Hp07_norm Data_to_plot_Hp10_norm Alldays_2020 C C2020 firstHp07_2020 ic ic2020 time1 time2 v v2020 number_of_days number_of_2020days 
     clear Siste*
     
     cd(global_raw_path)
end
% close(fig1)

clear global_N_File o ansatte_navn var10 var11 str_text 


 
%% 2. Nå sjekk de nye EPDs som %ender i *txt
%TODO: FOlder global_
% Excels with data location 
path=global_raw_path; 
cd(global_raw_path);
% With this there is no hidden files 
search_for_text=['*.txt'];
global_list_2=dir(search_for_text);

% Check in the folder you are 
d = dir('*.txt');% it shows hidden backup files also
N_File=numel(d)

% Delete this
% % Start active server for excel reading
% e = actxserver ('Excel.Application');
% h=waitbar(0,'Progress...');
%% Adjust to previous programs variables -todo matching soon
% global_path=path;
% global_path_figures=path_figures;
% global_B = datetime(2020,01,01);
%global_list_2=d2;

for o = 1:N_File %Loop start 
    qs=numel(global_list_2)-o;
    %clc;
    fprintf('Please o=wait %d second ',qs)
    fprintf(1,repmat('\n',1,1));
    
    %r = xlsread(d(o).name,2);
    cd(global_raw_path);
    excelfile_name=global_list_2(o).name;
   ansatte_navn=global_list_2(o).name(1:end-4);
    disp(['Working on: ' global_list_2(o).name])
    
    % Start active server for excel reading
    e = actxserver ('Excel.Application');
%     h=waitbar(0,'Progress...');
%     waitbar(o/N_File,h);
    ExcelWorkbook = e.workbooks.Open(fullfile(global_raw_path,global_list_2(o).name));

    %Sheet=ExcelWorkbook.Sheets.Item(2);
    
    exlFile=ExcelWorkbook;
    % We choose the first excel sheet
    exlSheet1 = exlFile.Sheets.Item(1);
    robj = exlSheet1.Columns.End(4);       % Find the end of the column
    numrows = robj.row;                    % And determine what row it is
    dat_range = ['A1:BA' num2str(numrows)]; % Read to the last row
    % cAN I HAVE IT AUTO? THE COLUMNS RANGE?
    rngObj = exlSheet1.Range(dat_range);

    % The contents of the sheet are here (our data)    
    exlData = rngObj.Value;

    %% Create table with correct headers
    C=exlData(1:end, :);
    T_temp=cell2table(C);
    T2=T_temp(:,1:5);
    
    %headersfromexcel={"Index"	"ProfileDateTime"	"TimeSuspect"	"EventString"	"Hp10"	"Hp
    T2.Properties.VariableNames={'index', 'ProfileDateTime', 'Hp10', 'Hp07', 'alarmtype'};
    cell_names=T2.Properties.VariableNames;
    % >super kul
    index_Hp07_name=find(strcmp(cell_names, 'Hp07'));
    index_Hp10_name=find(strcmp(cell_names, 'Hp10'));
    
    % Close the workbook now
    delete(ExcelWorkbook)
    
    clear T_temp
    SizeTable= size(T2);
    rows=SizeTable(1);
    
    for j=1:rows
        temp_time(:,j)=datetime(T2.ProfileDateTime(j), 'InputFormat', 'dd.MM.uuuu HH:mm:ss');
%         dataHp07(:,j)=cell2mat(T2.Hp07(j));
    end
    time_vec=temp_time';
    clear j 
    clear temp_time; 
    
    %% Cut out dates before 2020 and Non existing values in one shot. 
    idxHp07 = find(time_vec>global_B & T2.Hp07>0); %new
    idxHp10=find(time_vec>global_B & T2.Hp10>0); 

    Data_to_plot_Hp07(:,:) =(1/1000)*(T2.Hp07(idxHp07)); 
    time_vec_to_plot_Hp07(:,:) =time_vec(idxHp07);
    
    Data_to_plot_Hp10(:,:) = (1/1000)*T2.Hp10(idxHp10); 
    time_vec_to_plot_Hp10(:,:) =time_vec(idxHp10);
   
    %% FIGURE just to see if we got it right - and we did!!!
    %fig3=figure; plot(time_vec_to_plot_Hp10, Data_to_plot_Hp10, 'b*'); hold on; plot(time_vec_to_plot_Hp07, Data_to_plot_Hp07, 'r.');
    
   %% Normalize here start 
   %% Normalize extracting the initial value of 01.01.2020 (new 18012020)
   firstHp10_2020=T2.Hp10(idxHp10(1));
   firstHp07_2020=T2.Hp07(idxHp07(1));
   
   size_norm_Hp10=length(idxHp10);  %237
   size_norm_Hp07=length(idxHp07); % 176
    
   % make a normalized column Hp07
   for kk=1:size_norm_Hp07
       index_norm_temp=idxHp07(kk);
        T2.Hp07norm(index_norm_temp)=T2.Hp07(index_norm_temp)-firstHp07_2020;
   end
    clear kk index_norm_temp
     % make a normalized column Hp10
   for kk=1:size_norm_Hp10
       index_norm_temp=idxHp10(kk);
        T2.Hp10norm(index_norm_temp)=T2.Hp10(index_norm_temp)-firstHp10_2020;
   end
    clear kk index_norm  
    

  % Create Data_to_plot_norm variables 
  Data_to_plot_Hp10_norm(:,:)=T2.Hp10norm(idxHp10);
  Data_to_plot_Hp07_norm(:,:)=T2.Hp07norm(idxHp07);
  
  % Check that first value of normalized values of idx(1) is zero 
  
  %fig4=figure; plot(Data_to_plot_Hp10_norm, '*')
  
  % Now clear idx
  clear size_norm firstHp10_2020 first_Hp07_2020
  % End normalize
  
   
   % Normalize here end
   %% find total number of days above 01.01.2020
   % we do it only for Hp07
    v2020=datevec(time_vec(idxHp07));
    Alldays_2020=datetime(v2020(:,1:3));
    [C2020, ic2020]=unique(Alldays_2020);
    number_of_2020days=length(C2020);
    T2.total_days2020(1,1)=number_of_2020days;
    % end
    fprintf([global_list_2(o).name ':Total days in 2020:' num2str(number_of_2020days)]);
   
   
   %% save data to struct for SumDose
   % Siste data: 
  Siste_Hp07=Data_to_plot_Hp07_norm(end);
  Siste_Hp10=Data_to_plot_Hp10_norm(end);
  
  SumDoseMatlab.(ansatte_navn).SisteHp07=(1/1000)*Siste_Hp07;
  SumDoseMatlab.(ansatte_navn).SisteHp10=(1/1000)*Siste_Hp10;
  
  SumDoseMatlab.(ansatte_navn).totaldays=number_of_2020days;
   
    %% FIGURE 
    fig=figure; plot(time_vec_to_plot_Hp10, (1/1000)*Data_to_plot_Hp10_norm, 'b*'); hold on; plot(time_vec_to_plot_Hp07, (1/1000)*Data_to_plot_Hp07_norm, 'r.');
    
    %% For the figures
    % Make figures better
     legend('Hp10 norm','Hp07 norm','Location','northwest'); 
     xlabel('time');
     ylabel('Hp10(mSv) og Hp07(mSv)');
      
    % xlim([01-Jan-2020 05:00:00 27-Nov-2020 20:00:00])
     time1=datetime("01.01.2020 00:00:09",'InputFormat', "dd.MM.yyyy HH:mm:ss");
     time2=datetime("27.11.2020 14:00:00", 'InputFormat', "dd.MM.yyyy HH:mm:ss");
     xlim([time1 time2]); 
     ylim([0 6])
     %% For anonymt plot
     %title(['Hp10 ansatt nr '  num2str(o)]);
     
     %% For ikke anonymt plot
     title([ansatte_navn]);
     
     %% Save figures
     savename=[excelfile_name(1:end-5) '_results'];
     cd(global_path);
     cd(global_figures_path)
     mkdir(savename)
     cd(savename)
     savefig(fig,savename)
     saveas(fig, savename, 'tiff')
     saveas(fig, savename, 'fig')
     
     %% Add figure to rapport 
     fig1=Figure(fig);
     fig1.Width="4.2in";
     fig1.Height="4in";
     %fig1.Scaling="custom";
    
     
       %% Create chapters
    ch=Chapter;
    ch.Title=ansatte_navn;
    ch.Numbered=true; 
   % ch.number=16+o;
    add(rpt, ch);
     add(rpt, fig1); 
     close(fig)
     
     
     cd(global_raw_path);
     save T2 T2
     
     clear SizeTable rows Data_to_plot_Hp07 Data_to_plot_Hp10 idx time_vec time_vec_to_plot T2 excelfile_name savename qs rows columns
     clear time_vec_to_plot_Hp10 time_vec_to_plot_Hp07 Data_to_plot_Hp10 Data_to_plot_Hp10_norm Data_to_plot_Hp07_norm index_Hp07_name index_Hp10_name indices_Hp07 indices_Hp10
    clear idxHp07 idxHp10 size_norm_Hp07 size_norm_Hp10 index_norm_temp
     cd(global_path)
     % Close Excel 
     %quit(e)
     clear e
     

end

ch=Chapter;
ch.Title="Agnell Kristina";
ch.Numbered=true; 
add(rpt, ch);
T=Text("Det var ikke noe data for Agnell Kristina");
T.FontSize="12pt";
T.Italic=true;
add(rpt, T);

%% Final end; nå kan du lukke report
close(rpt);

save SumDoseMatlab SumDoseMatlab

xticksmaria=cell(17,1);
global_list_all=[global_list; global_list_2];
for lole=1:length(global_list_all), 
    xticksmine3{lole,1}=global_list_all(lole).name(1:end-4); 
end

fig2=figure;
for kkj=1:length(global_list_all)
hold on;
ansatte_navn=global_list_all(kkj).name(1:end-4);
bar(kkj, SumDoseMatlab.(ansatte_navn).SisteHp07, 'g');
hold on; bar(kkj,SumDoseMatlab.(ansatte_navn).SisteHp10, 'y');
%hold on; bar(kkj, SumDoseMatlab.(ansatte_navn).totaldays, 'm');
end
hold off
ax = gca;
ax.XTick = [1:1:length(global_list_all)];
% xticksmaria={[xticksmaria; global_list_all(kkj).name(1:end-4)]}
ax.XTickLabels = xticksmine3;
ax.XTickLabelRotation = 45;

% Make figures better
     legend('Hp10 norm','Hp07 norm','Location','northwest'); 
     xlabel('Ansatt');
     ylabel('Hp10(mSv) og Hp07(mSv)');

savename=['SumDose_results'];
     cd(global_path);
     cd(global_figures_path)
     mkdir(savename)
     cd(savename)
     savefig(fig2,savename)
     saveas(fig2, savename, 'tiff')
     saveas(fig2, savename, 'fig')
     
    