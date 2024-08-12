%__________________________________________________________________________
%
% REEFMOD-GBR combined Hindcast (2008-2023.5) and Forecast (2024-2100)
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 12/2019. Last update: 02/2024
%__________________________________________________________________________


%% INITIAL REEF STATES WITH AIMS LTMP TRANSECT AND MANTA TOW DATA
load('AIMS_LTMP_Transect_means.mat') 
% Complete (1992-2017) LTMP transect dataset at 5-6m averaged per reef (2-3 sites per reef).
% Total cover (TOT) + cover of each coral group (G1....G6).
% Reefs assigned to Karlo's polygons based on shortest distance between survey and centroid coordinates
% (see 'RE-ARRANGE_LTMP_TRANSECTS.m'). There are discrepancies in reef names and shelf position 
% between AIMS database and Karlo's list of 3806. using Karlo's designations here for consistency
% 'YEAR_Reefmod' is modelled year (surveys assigned to each season from their actual dates) 

%% 1) Build the average community compo for the corresponding years
survey_year = LTMP_Transects_means.YEAR_Reefmod;
select_LTMP_compo = find(survey_year>=2004 & survey_year<2009); % select transect data from 2004 to 2008 to initialise community compo

TOT_COVER = table2array(LTMP_Transects_means(select_LTMP_compo,9)); % 341 reefs surveyed between 2004-08 (123 unique reefs)
Z = table2array(LTMP_Transects_means(select_LTMP_compo,10:15))./TOT_COVER(:,ones(6,1)); % relative proportion of each group (each row sums to 1)
LTMP_Transects_compos = [LTMP_Transects_means(select_LTMP_compo,[16 17 3]) array2table(Z)]; % Extract AIMS sector (3), not Karlo's (18) as they differ
LTMP_Transects_compos.NEW_SHELF(LTMP_Transects_compos.SHELF=='I')=1;
LTMP_Transects_compos.NEW_SHELF(LTMP_Transects_compos.SHELF=='M')=2;
LTMP_Transects_compos.NEW_SHELF(LTMP_Transects_compos.SHELF=='O')=3;

% AIMS shelf designation needs the following corrections to match ReefMod shelf position
LTMP_Transects_compos(LTMP_Transects_compos.KarloID==684,10)=array2table(2);
LTMP_Transects_compos(LTMP_Transects_compos.KarloID==836,10)=array2table(2);
LTMP_Transects_compos(LTMP_Transects_compos.KarloID==1034,10)=array2table(2);
LTMP_Transects_compos(LTMP_Transects_compos.KarloID==1051,10)=array2table(2);
LTMP_Transects_compos(LTMP_Transects_compos.KarloID==1885,10)=array2table(2);
LTMP_Transects_compos(LTMP_Transects_compos.KarloID==1887,10)=array2table(2);
LTMP_Transects_compos(LTMP_Transects_compos.KarloID==1911,10)=array2table(2);
LTMP_Transects_compos(LTMP_Transects_compos.KarloID==2230,10)=array2table(1);
LTMP_Transects_compos(LTMP_Transects_compos.KarloID==3321,10)=array2table(3);
LTMP_Transects_compos(LTMP_Transects_compos.KarloID==3629,10)=array2table(2);
LTMP_Transects_compos(LTMP_Transects_compos.KarloID==3630,10)=array2table(2);

omean = @(x) mean(x,'omitnan');
% Calculate average community composition (relative props) for each sector(1-11)*shelf positon(1-3)
AIMS_AVERAGECOMPO = table2array(varfun(omean, LTMP_Transects_compos,'GroupingVariables',{'SectorCode','NEW_SHELF'},...
    'InputVariables',{'Z1','Z2','Z3','Z4','Z5','Z6'})); % 20 combinations of sector*shelf (out of 33)

% Add the missing sector/shelf combinations
% For missing compo take the one of the preceding sector (north) for the same shelf position
SECTOR = AIMS_AVERAGECOMPO(:,1);
SHELF = AIMS_AVERAGECOMPO(:,2);

% Note there is no inshore reef in the Swain sector!
Missings = [1 1 NaN AIMS_AVERAGECOMPO(SECTOR==4&SHELF==1,4:9) ;...
            1 2 NaN AIMS_AVERAGECOMPO(SECTOR==3&SHELF==2,4:9) ;...
            1 3 NaN AIMS_AVERAGECOMPO(SECTOR==3&SHELF==3,4:9) ;...
            2 1 NaN AIMS_AVERAGECOMPO(SECTOR==4&SHELF==1,4:9) ;...
            2 2 NaN AIMS_AVERAGECOMPO(SECTOR==3&SHELF==2,4:9) ;...
            2 3 NaN AIMS_AVERAGECOMPO(SECTOR==3&SHELF==3,4:9) ;...
            3 1 NaN AIMS_AVERAGECOMPO(SECTOR==4&SHELF==1,4:9) ;... 
            5 3 NaN AIMS_AVERAGECOMPO(SECTOR==4&SHELF==3,4:9) ;...          
            7 1 NaN AIMS_AVERAGECOMPO(SECTOR==6&SHELF==1,4:9) ;...
            7 2 NaN AIMS_AVERAGECOMPO(SECTOR==6&SHELF==2,4:9) ;...
            7 3 NaN AIMS_AVERAGECOMPO(SECTOR==6&SHELF==3,4:9) ;...
            9 3 NaN AIMS_AVERAGECOMPO(SECTOR==8&SHELF==3,4:9) ];
            
AIMS_AVERAGECOMPO = array2table([AIMS_AVERAGECOMPO ; Missings]) ; % Now 32 -> all possible sector*shelf combinations populated

% Assign to each of the 3806 reefs the average community compo representative of shelf and sector
relative_cover = table2array(outerjoin(GBR_REEFS,AIMS_AVERAGECOMPO, 'LeftKeys', [10 11],'RightKeys', [2 1], 'LeftVariables',[1 10 11], 'RightVariables',[4:9]));

[~,I]=sort(relative_cover(:,1));
relative_cover = relative_cover(I,:); % sort by increasing reef ID

% Then swap for the surveyed reefs (n=123) which compo was actually recorded (averaged over the selected period (2006-2008)
AIMS_AVERAGECOMPO_REEF = table2array(varfun(omean, LTMP_Transects_compos,'GroupingVariables',{'KarloID'},...
    'InputVariables',{'Z1','Z2','Z3','Z4','Z5','Z6'}));

relative_cover(AIMS_AVERAGECOMPO_REEF(:,1),4:end) = AIMS_AVERAGECOMPO_REEF(:,3:end); % including 123 reefs with actual community composition recorded

%% 2) Assign to each of the 3806 reefs the average total cover representative of shelf and sector between 2006-2008
StartDate = 2006 ; % earliest date used to initialize coral cover
EndDate = 2009 ; % exclusive

% Calculate average total cover of reefs surveyed with transects %%%%%%%%%%%%%%%%%%%
select_TR = find(LTMP_Transects_means.YEAR_Reefmod >= StartDate & LTMP_Transects_means.YEAR_Reefmod < EndDate);
AIMS_TR = LTMP_Transects_means(select_TR,[16 17 3 9 5]); %5: Years, 10: cover tot, 9: YEAR_reefmod
AIMS_TR.Properties.VariableNames = {'Reef_ID','SECTOR','AIMS_SHELF','CCOVER','YEAR'};

AIMS_TR.SHELF(AIMS_TR.AIMS_SHELF=='I')=1;
AIMS_TR.SHELF(AIMS_TR.AIMS_SHELF=='M')=2;
AIMS_TR.SHELF(AIMS_TR.AIMS_SHELF=='O')=3;
AIMS_TR = AIMS_TR(:,[1 2 6 4 5]);
% length(unique(AIMS_TR.Reef_ID)); %123 reefs surveyed using transects between 2006-2008

% AIMS shelf designation needs the following corrections to match ReefMod shelf position
AIMS_TR(AIMS_TR.Reef_ID==684,3)=array2table(2);
AIMS_TR(AIMS_TR.Reef_ID==836,3)=array2table(2);
AIMS_TR(AIMS_TR.Reef_ID==1034,3)=array2table(2);
AIMS_TR(AIMS_TR.Reef_ID==1051,3)=array2table(2);
AIMS_TR(AIMS_TR.Reef_ID==1885,3)=array2table(2);
AIMS_TR(AIMS_TR.Reef_ID==1887,3)=array2table(2);
AIMS_TR(AIMS_TR.Reef_ID==1911,3)=array2table(2);
AIMS_TR(AIMS_TR.Reef_ID==2230,3)=array2table(1);
AIMS_TR(AIMS_TR.Reef_ID==3321,3)=array2table(3);
AIMS_TR(AIMS_TR.Reef_ID==3629,3)=array2table(2);
AIMS_TR(AIMS_TR.Reef_ID==3630,3)=array2table(2);

% Select manta tow surveys %%%%%%%%%%%%%%%
% Feb 2024: now using AIMS + FMP manta tows
load('ALL_MantaTow_CORAL_1991_2024.mat')
select_MT = find(MantaTow_CORAL_1991_2024.YearReefMod >= StartDate & MantaTow_CORAL_1991_2024.YearReefMod < EndDate); % 258 reef surveys
ALL_MT = MantaTow_CORAL_1991_2024(select_MT,[7 12 5 6]); % Note this is using ReefMod shelf position (different than AIMS shelf classification)
ALL_MT.Properties.VariableNames = {'Reef_ID','SHELF','MEAN_LIVE_CORAL','YEAR'}; % Note Reef_ID not available here

% convert Manta tows into transect equivalent (calibration model updated Jan 2024) %%%%%%%%%%%%%
load('LTMP_Transect2Tow_2024.mat') % linear model allowing to predict transect coral (y) cover from manta-tow coral cover (x)
% Model in the form sqrt(y)~sqrt(x) (ie, with forced 0 intercept)
ALL_MT.MEAN_LIVE_CORAL_EQ = predict(LTMP_Transect2Tow_Model2,sqrt(ALL_MT.MEAN_LIVE_CORAL)).^2;
ALL_MT.SECTOR = GBR_REEFS.AIMS_sector(ALL_MT.Reef_ID);
ALL_MT = ALL_MT(:,[1 6 2 5 4]);
ALL_MT.Properties.VariableNames = {'Reef_ID','SECTOR','SHELF','CCOVER','YEAR'};

% Bind transect and tow data
ALL_CORALtmp = [ALL_MT ; AIMS_TR];
ALL_CORAL = varfun(omean, ALL_CORALtmp,'GroupingVariables',...
    {'Reef_ID';'SECTOR';'SHELF'}, 'InputVariables','CCOVER'); %average across sites/method

CORAL_AVERAGESTATE = table2array(varfun(omean, ALL_CORAL,'GroupingVariables',...
    {'SECTOR','SHELF'}, 'InputVariables',{'Fun_CCOVER'})); % 25 combinations of sector*shelf (out of 32)
% Exclude sectors represented by less than 3 reefs
CORAL_AVERAGESTATE = CORAL_AVERAGESTATE(CORAL_AVERAGESTATE(:,3)>=3,:); % 19 combinations of sector*shelf (out of 32)

% Add the missing sector/shelf combinations
SECTOR = CORAL_AVERAGESTATE(:,1);
SHELF = CORAL_AVERAGESTATE(:,2);

% Note there is no inshore reef in the Swain sector!
Missings = [1 1 NaN CORAL_AVERAGESTATE(SECTOR==4&SHELF==1,4) ;...
    1 2 NaN CORAL_AVERAGESTATE(SECTOR==3&SHELF==2,4) ;...
    1 3 NaN CORAL_AVERAGESTATE(SECTOR==3&SHELF==3,4) ;...
    2 1 NaN CORAL_AVERAGESTATE(SECTOR==4&SHELF==1,4) ;...
    2 2 NaN CORAL_AVERAGESTATE(SECTOR==3&SHELF==2,4) ;...
    2 3 NaN CORAL_AVERAGESTATE(SECTOR==3&SHELF==3,4) ;...
    3 1 NaN CORAL_AVERAGESTATE(SECTOR==4&SHELF==1,4) ;...  
    5 3 NaN CORAL_AVERAGESTATE(SECTOR==4&SHELF==3,4) ;...
    7 1 NaN CORAL_AVERAGESTATE(SECTOR==6&SHELF==1,4) ;...
    7 3 NaN CORAL_AVERAGESTATE(SECTOR==6&SHELF==3,4) ;...
    9 1 NaN CORAL_AVERAGESTATE(SECTOR==8&SHELF==1,4) ;...
    9 3 NaN CORAL_AVERAGESTATE(SECTOR==8&SHELF==3,4) ;...
    11 2 NaN CORAL_AVERAGESTATE(SECTOR==9&SHELF==2,4) ] ;
        
CORAL_AVERAGESTATE = array2table([CORAL_AVERAGESTATE ; Missings],'VariableNames', {'Sector';'Shelf';'N';'TOT'} ) ;

init_coral = table2array(outerjoin(GBR_REEFS,CORAL_AVERAGESTATE, 'LeftKeys', [10 11],'RightKeys', [2 1], 'LeftVariables',[1 10 11], 'RightVariables',[4]));
init_coral = init_coral(I,:); % sort by increasing reef ID (just in case)

clear Z TOT_COVER survey_year select_nan omean SECTOR SHELF Missings ALL_CORALtmp

% Then for each surveyed reef swap the assigned average by the actual value of total cover
% This gives x reefs initialised with the total cover observed, and y reefs initialised with a regional average
init_coral(ALL_CORAL.Reef_ID,4) = ALL_CORAL.Fun_CCOVER;


%% 3) Now populate each reef with a random cover (normally distributed around specified mean)
X = init_coral(META.reef_ID,4);
Y = relative_cover(META.reef_ID,4:end);
varSDcoral = 0.2*ones(1,META.nb_coral_types);

init_coral_cover = X(:,ones(size(Y,2),1)).* Y ;

for s=1:6
    init_coral_cover(:,s) = normrnd(init_coral_cover(:,s), varSDcoral(1,s)*init_coral_cover(:,s))/100;
end

init_coral_cover(init_coral_cover<0.005)=0.005; % minimum 0.5% cover for each group
TCmax = 0.8; %maximum initial total cover is 80%
select_TCmax = find(sum(init_coral_cover,2)>TCmax);
adjust_TCmax = sum(init_coral_cover(select_TCmax,:),2);
init_coral_cover(select_TCmax,:)=TCmax*init_coral_cover(select_TCmax,:)./adjust_TCmax(:,ones(size(init_coral_cover,2),1));

%% Now populate random rubble and sand covers
init_rubble = 0.1; % 10% on average 
init_sand = 0.3; % 30% sand on average

varSDother = 0.2;

init_rubble_cover = normrnd(init_rubble*100, varSDother*init_rubble*100, length(META.reef_ID), 1)/100 ;
init_rubble_cover(init_rubble_cover<0.01) = 0.01;

% Update Nov 2022: now sand determined by Roelfsema et al (2021)'s geomorphic and benthic maps
init_sand_cover = normrnd(GBR_REEFS.UNGRAZABLE(META.reef_ID)*100, varSDother*GBR_REEFS.UNGRAZABLE(META.reef_ID)*100)/100 ;
init_sand_cover(init_sand_cover<0.05) = 0.05; % impose minimum of 5%

% Adjust if too high (coral + sand > 95%). We leave 5% free space for a safe intialisation
CHECK = sum(init_coral_cover,2) + init_sand_cover;
init_sand_cover(CHECK>0.95) = 0.95 - sum(init_coral_cover(CHECK>0.95,:),2);

clear init_coral init_sand init_rubble relative_cover X Y varSD adjust_TCmax select_TCmax TCmax varSDother CHECK
clear AIMS_AVERAGECOMPO AIMS_AVERAGECOMPO_REEF CORAL_AVERAGESTATE AIMS_TR LTMP_Transects_compos LTMP_Transects_means
clear select_LTMP_compo select_MT select_TR varSDcoral MantaTow_CORAL_1991_2024 LTMP_Transect2Tow_Model2 ALL_CORAL ALL_MT

%% COTS densities - INPUT DENSITIES MUST BE PER GRID (400m2), NOT PER TOW
REEF_COTS.densities_M = nan(META.nb_reefs,META.nb_time_steps+1);
REEF_COTS.densities_SD = nan(META.nb_reefs,META.nb_time_steps+1);

% First, initialize CoTS populations with CoCoNet hindcast densities of adult CoTS (mean and SD per grid) provided by Scott
load('GBR_past_COTS_CoCoNet.mat')
% Matrices of mean density and SD already converted into nb CoTS per grid for the period 1985-2017.
% Values are mean based on 50 replicate runs, interpolated from 3638 reefs - see 'REEF_COTS_INTERPOLATION_1985_2017.m'
% So 33 years * 2 seasons = 66 columns
CoCoNet_years = 1985:0.5:2017.5;

start_year = find(CoCoNet_years==2007.5) ; %first column is 1985, so 2007 is 1+22
REEF_COTS.densities_M(:,1) = PAST_COTS_M_COCONET(META.reef_ID,start_year); % data are already converted into CoTS densities per grid
REEF_COTS.densities_SD(:,1) = PAST_COTS_SD_COCONET(META.reef_ID,start_year); % data are already converted into CoTS densities per grid

% Then we force CoTS density to past observations where available
% load('GBR_PAST_COTS_1992_2020.mat') % Mean coTS per tow between 1992 and 2000 (AIMS LTMP and FMP)
load('GBR_PAST_COTS_1992_2023.mat') % Mean coTS per tow between 1992 and 2023 (AIMS LTMP, FMP and Control Program)
% This is the average number of CoTS per manta tow. Missing values are now identified as NaN, not -1 as previously.

GBR_PAST_COTS_PER_GRID = (GBR_PAST_COTS_NEW(META.reef_ID,:)/0.22)*(1500/2500); % Conversion: 1500 COTS per km2 ~0.22 per manta tow (Moran and De'ath 1986)
% which is further divided by 2500 to get number of COTS per 400m2 (reef grid size)

COTS_years = 1992:0.5:2023;
start_step = find(COTS_years==2007.5);
end_step = size(GBR_PAST_COTS_PER_GRID,2);
length_history = length(start_step:end_step);

RECENT_COTS = max(GBR_PAST_COTS_PER_GRID(:,(start_step-3):start_step),[],2); % maximum density observed from 2006 to 2007.5
tmp2 = find(isnan(RECENT_COTS)==0);

% Replace CoCoNet predictions with recent observations (2006-2007.5)
REEF_COTS.densities_M(find(isnan(RECENT_COTS)==0),1)= RECENT_COTS(find(isnan(RECENT_COTS)==0));
REEF_COTS.densities_SD(find(isnan(RECENT_COTS)==0),1)= 0.28 * RECENT_COTS(find(isnan(RECENT_COTS)==0)); 
% (with 0.28 = slope of linear model between SD and mean CoCoNet predictions)

% Populate from 2008 onwards with observations
REEF_COTS.densities_M(:,2:length_history) = GBR_PAST_COTS_PER_GRID(:,(start_step+1):end_step) ;
REEF_COTS.densities_SD(:,2:length_history) = 0.28 * GBR_PAST_COTS_PER_GRID(:,(start_step+1):end_step) ;

clear GBR_PAST_COTS about_GBR_COTS PAST_COTS_M_COCONET PAST_COTS_SD_COCONET CoCoNet_years end_step
clear GBR_PAST_COTS_NEW GBR_PAST_COTS_NEW_ID GBR_PAST_COTS_PER_GRID about_GBR_PAST_COTS MantaTow_1985_2020

%% SCENARIO OF THERMAL STRESS
if META.doing_genetics == 0 % NO GENETIC ADAPTATION (WITH THERMAL OPTIMUM)
    
    DHW = zeros(length(META.reef_ID),META.nb_time_steps);
    
    %% PAST THERMAL STRESS REGIME
    % Using NOAA Coral Reef Watch 5km product: max DHW every year from 1985 to 2023 from closest 5x5 km pixel (updated Jan 2024)   
    load('GBR_past_DHW_CRW_5km_1985_2023.mat')
    GBR_PAST_DHW = GBR_PAST_DHW(:,24:39); %select 2008 to 2023
  
    end_hindcast = size(GBR_PAST_DHW,2);
    DHW(:,1:2:end_hindcast*2) = GBR_PAST_DHW(META.reef_ID,:);
    clear GBR_PAST_DHW

    %% FUTURE THERMAL STRESS (CMIP6)
    start_future = 33; % Forecast starts in 2024 (step 33)
    % Need to be odd number (always start future in summer)
    
    if META.nb_time_steps >= start_future
        
        % The DHW matrices start in 2000 (column 6) so year 2024 is column 30
        
        if OPTIONS.SSP=='119' && ismember(OPTIONS.GCM, {'GFDL-ESM4', 'MPI-ESM1-2-HR', 'NorESM2-LM'})==1           
            error('##### REEFMOD ERROR #### SSP1-1.9 is not available for the climate models FDL-ESM4, MPI-ESM1-2-HR and NorESM2-LM')     
        end
        
        % Load the selected forecast scenario of DHW
        load([char(OPTIONS.GCM) '_' char(OPTIONS.SSP) '_annual_DHW_max.mat'])
       
        if META.nb_time_steps-start_future-1 > 2*size(max_annual_DHW(:,30:end),2)   
            error('##### REEFMOD ERROR #### Simulated timeframe inconsistent with the DHW forecast. NB_TIME_STEPS has to be <= 31+155')
        else
            
            % Select the specified timeframe in the available forecast
            DHW_FORECAST = max_annual_DHW(META.reef_ID,30:(30+(META.nb_time_steps-start_future-1)/2));
            % Let's shuffle available years within each decade
            DHW_FORECAST_shuffled = nan(size(DHW_FORECAST));
            start = 1; % first year of the selected forecast
            remain = size(DHW_FORECAST,2); % number of years still available
            
            while remain > 10
                
                sample = randperm(10); % sample at random within the decade
                DHW_FORECAST_shuffled(:,start:(start-1+length(sample)))= DHW_FORECAST(:,start-1+sample); % assign the shuffled years
                start = start + length(sample);
                remain = remain - length(sample);
            end
            
            % Last years available to be shuffled as well (Less than a decade s available)
            sample = randperm(remain);
            DHW_FORECAST_shuffled(:,start:(start-1+length(sample)))= DHW_FORECAST(:,start-1+sample);
            
            % Finally assign to the DHW matrix (only in summers)
            DHW(:,start_future:2:end) = DHW_FORECAST_shuffled;
        end
    end
    
else
    error('Cannot run genetic adaptation with this version of ReefMod')
end


%% PAST STORM REGIME 2008 to 2023
% Category cyclones (Saffir-Simpson scale) assigned to the spatial prediction of >4m wave height from Puotinen et al (2016).
% Original set of predictions is matrice of 0/1 compiled by Rob for the entire GBR (3806 reefs). 
% Then, occurence of damaging wave is blended with category cyclone estimated from maximum sustained wind and distance to cyclone track
% (with 10% increase to meet US standards) measured along each real cyclone track (data extracted from BoM Database of past cyclone tracks).
load('GBR_cyclones_2008-2020_NEW.mat')

CYCLONE_CAT = zeros(length(META.reef_ID),2*size(GBR_PAST_CYCLONES_NEW,2)) ;
CYCLONE_CAT(:,1:2:end) = GBR_PAST_CYCLONES_NEW(META.reef_ID,:); % 26 columns from 2008 to 2020 (summer/winter)
clear GBR_PAST_CYCLONES_NEW

%% Update Jan 2024: PATCH FOR 2021, 2022 & 2023 (assuming no severe cyclone)
CYCLONE_CAT = [CYCLONE_CAT zeros(length(META.reef_ID), 2*3)];
    
%% FUTURE STORM REGIME (from 2023 onwards)
if META.nb_time_steps > 31
    
    load('Reef_Cyclone_TimeSeries_Count_Cat.mat')
    
    FUTURE_CYCLONE_CAT =  zeros(length(META.reef_ID),2*size(Cyc_cat,2)) ;
    FUTURE_CYCLONE_CAT(:,1:2:end) = Cyc_cat(GBR_REEFS.Nick_ID(META.reef_ID),:,simul);
    
    CYCLONE_CAT = [CYCLONE_CAT FUTURE_CYCLONE_CAT];
    
end

clear GBR_DHM_Max DHM_Max FUTURE_COTS_M FUTURE_COTS_SD CURRENT_CORAL_COVER CURRENT_COTS_DENSITIES
clear CURRENT_RUBBLE_COVER CURRENT_NONGRAZABLE ALL_CYCLONE_CAT Cyc_cat DHWtmp Model_reef_indices ModelReef_ID ModelReef_ID_label
clear SST_GBRtmp Years_unique yr nb_hindcast_runs rel T_opt
clear nb_time_steps Reef_SST_Mean_Yr total_coral_cover Yr Topt_scenario_YM varSD resize_Topt_list
clear FUTURE_CYCLONE_CAT AIMS_ALL AIMS_ALLtmp AIMS_MT AIMS_TR max_annual_DHW LTMP_Transect2Tow_Model
