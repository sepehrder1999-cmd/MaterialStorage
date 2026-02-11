function model=CreateModel()

    % Table1: Geometry and Time Data of Site Facilities
    % Site Data
    site_width = 55;
    site_height = 75;

    model.site_width=site_width;
    model.site_height=site_height;

    
    Stage=[1,2,3];

    % Duration in weeks
    StageDuration=[18,26,26];
    st=numel(Stage);
    
    for x=1:st
        model.Stage(x).Name = ['Stage ' num2str(x)];
        model.Stage(x).Duration = StageDuration(x);
    end
    
    % Fixed Facilities
    % Building B1 dimensions and position
    B1.x = 5;
    B1.y = 45;
    B1.Lx = 45;
    B1.Ly = 25;
    B1.activestages = [1, 2, 3];
    B1.Fixed = 1;
    Facilities.B1=B1;
    
    % Building B2 dimensions and position
    B2.x = 30;
    B2.y = 15;
    B2.Lx = 20;
    B2.Ly = 20;
    B2.activestages = [2, 3];
    B2.Fixed=1;
    Facilities.B2=B2;

    % Access Gate dimensions and position
    G.name="Gate";
    G.x = 0;
    G.y = 0;
    G.Lx = 10;
    G.Ly = 1;
    G.activestages = [1, 2, 3];
    G.Fixed=1;
    Facilities.G=G;

    % Tower crane
    F1.name = "Tower crane";
    F1.Lx=8;
    F1.Ly=8;
    F1.activestages = [1, 2, 3];
    F1.Type='S';
    F1.RelCost=0;
    Facilities.F1=F1;

    % Site office trailor (1)
    F2.name="Site office trailor (1)";
    F2.Lx=14;
    F2.Ly=4;
    F2.activestages = [1, 2, 3];
    F2.Type='M';
    F2.RelCost=6000;
    Facilities.F2=F2;

     % Site office trailor (2)
    F3.name = "Site office trailor (2)";
    F3.Lx=11;
    F3.Ly=3;
    F3.activestages = [2, 3];
    F3.Type='M';
    F3.RelCost=4000;
    Facilities.F3=F3;

    %  Fabrication Area
    F4.name = "Fabrication Area";
    F4.Lx=15;
    F4.Ly=10;
    F4.activestages = [1, 2, 3];
    F4.Type='M';
    F4.RelCost=2000;
    Facilities.F4=F4;

    % Dump Area
    F5.name = "Dump Area";
    F5.Lx=15;
    F5.Ly=15;
    F5.activestages = [1, 2];
    F5.Type='M';
    F5.RelCost=0;
    Facilities.F5=F5;

    % Lay-down Area
    F6.name = "Lay-down Area";
    F6.Lx=10;
    F6.Ly=10;
    F6.activestages = [1, 2, 3];
    F6.Type='M';
    F6.RelCost=3000;
    Facilities.F6=F6;

    % Labor Rest Area
    F7.name = "labor rest Area";
    F7.Lx=5;
    F7.Ly=5;
    F7.activestages = [1, 2, 3];
    F7.Type='M';
    F7.RelCost=500;
    Facilities.F7=F7;

    % Fuel Storage
    F8.name = "labor rest Area";
    F8.Lx=3;
    F8.Ly=3;
    F8.activestages = [1, 2, 3];
    F8.Type='M';
    F8.RelCost=100;
    Facilities.F8=F8;


    model.Facilities = Facilities; 

    % Travel Cost Rates ($/m) Table 2

    % Define the list of facilities in order
    FacilityList = [fieldnames(model.Facilities)];

    n = length(FacilityList);

    % Create a map from facility name to index
    facilityMap = containers.Map(FacilityList, 1:n);

    % Initialize cost matrix with NaN (undefined connections)
    travelCostRate = NaN(n, n);
    model.facilityIndex = {'B1','B2','G','F1','F2','F3','F4','F5','F6','F7'};


    % B1 row
    travelCostRate(facilityMap('B1'), facilityMap('B1')) = 0;
    travelCostRate(facilityMap('B1'), facilityMap('B2')) = 0;
    travelCostRate(facilityMap('B1'), facilityMap('G')) = 0;

    travelCostRate(facilityMap('B1'), facilityMap('F1')) = 150;
    travelCostRate(facilityMap('B1'), facilityMap('F2')) = 50;
    travelCostRate(facilityMap('B1'), facilityMap('F3')) = 50;
    travelCostRate(facilityMap('B1'), facilityMap('F4')) = 90;
    travelCostRate(facilityMap('B1'), facilityMap('F5')) = 20;
    travelCostRate(facilityMap('B1'), facilityMap('F6')) = 70;
    travelCostRate(facilityMap('B1'), facilityMap('F7')) = 15;

    % B2 row
    travelCostRate(facilityMap('B2'), facilityMap('B2')) = 0;
    travelCostRate(facilityMap('B2'), facilityMap('G')) = 0;

    travelCostRate(facilityMap('B2'), facilityMap('F1')) = 100;
    travelCostRate(facilityMap('B2'), facilityMap('F2')) = 40;
    travelCostRate(facilityMap('B2'), facilityMap('F3')) = 40;
    travelCostRate(facilityMap('B2'), facilityMap('F4')) = 60;
    travelCostRate(facilityMap('B2'), facilityMap('F5')) = 15;
    travelCostRate(facilityMap('B2'), facilityMap('F6')) = 40;
    travelCostRate(facilityMap('B2'), facilityMap('F7')) = 15;

    % G row (site gate)
    travelCostRate(facilityMap('G'), facilityMap('F1')) = 0;
    travelCostRate(facilityMap('G'), facilityMap('F2')) = 2;
    travelCostRate(facilityMap('G'), facilityMap('F3')) = 2;
    travelCostRate(facilityMap('G'), facilityMap('F4')) = 1;
    travelCostRate(facilityMap('G'), facilityMap('F5')) = 30;
    travelCostRate(facilityMap('G'), facilityMap('F6')) = 0;
    travelCostRate(facilityMap('G'), facilityMap('F7')) = 0;

    % F1 row
    travelCostRate(facilityMap('F1'), facilityMap('F2')) = 0;
    travelCostRate(facilityMap('F1'), facilityMap('F3')) = 0;
    travelCostRate(facilityMap('F1'), facilityMap('F4')) = 30;
    travelCostRate(facilityMap('F1'), facilityMap('F5')) = 4;
    travelCostRate(facilityMap('F1'), facilityMap('F6')) = 25;
    travelCostRate(facilityMap('F1'), facilityMap('F7')) = 0;

    % F2 row
    travelCostRate(facilityMap('F2'), facilityMap('F3')) = 20;
    travelCostRate(facilityMap('F2'), facilityMap('F4')) = 5;
    travelCostRate(facilityMap('F2'), facilityMap('F5')) = 0;
    travelCostRate(facilityMap('F2'), facilityMap('F6')) = 5;
    travelCostRate(facilityMap('F2'), facilityMap('F7')) = 0;

    % F3 row
    travelCostRate(facilityMap('F3'), facilityMap('F4')) = 5;
    travelCostRate(facilityMap('F3'), facilityMap('F5')) = 0;
    travelCostRate(facilityMap('F3'), facilityMap('F6')) = 5;
    travelCostRate(facilityMap('F3'), facilityMap('F7')) = 0;


    % F4 row
    travelCostRate(facilityMap('F4'), facilityMap('F5')) = 0;
    travelCostRate(facilityMap('F4'), facilityMap('F6')) = 30;
    travelCostRate(facilityMap('F4'), facilityMap('F7')) = 0;

    % F5 row
    travelCostRate(facilityMap('F5'), facilityMap('F6')) = 0;
    travelCostRate(facilityMap('F5'), facilityMap('F7')) = 0;

    % F6 row
    travelCostRate(facilityMap('F6'), facilityMap('F7')) = 0;

        % Set diagonal elements to 0 (cost from a facility to itself)
    for a = 1:n
        travelCostRate(a, a) = 0;
    end
        for a = 1:n
            for b = a+1:n
                travelCostRate(b,a) = travelCostRate(a,b);
            end
        end
    model.travelCostRate=travelCostRate;


    %  file that contains material demand flow over project timeline
    demandTable = readtable('MaterialDemand.xlsx');

    % Table 3 Storage Footprints of Construction Materials

    % Material M1: Rebar

    materials(1).ID = 'M1';
    materials(1).Name = 'Rebar';
    materials(1).Unit = 'Ton';

    materials(1).PurchaseCost = [
        struct('Stage', 1, 'Quantity', [0, 100], 'PCR', 650)
        struct('Stage', 1, 'Quantity', [100, 200], 'PCR', 550)
        struct('Stage', [2, 3], 'Quantity', [0, 100], 'PCR', 750)
        struct('Stage', [2, 3], 'Quantity', [100, 200], 'PCR', 650)
        ];

    materials(1).DeliveryCost = [
        struct('Quantity', [0, 25], 'DLC', 600)
        struct('Quantity', [25, 50], 'DLC', 1200)
        struct('Quantity', [50, 75], 'DLC', 1800)
        struct('Quantity', [75, 100], 'DLC', 2400)
        struct('Quantity', [100, 125], 'DLC', 3000)
        struct('Quantity', [125, 150], 'DLC', 3600)
        struct('Quantity', [150, 175], 'DLC', 4200)
        struct('Quantity', [175, 200], 'DLC', 4800)
        ];

    materials(1).StorageFootprint = [
        struct('Qrange', [0, 32], 'Lx', 15, 'Ly', 2)
        struct('Qrange', [32, 64], 'Lx', 15, 'Ly', 4)
        struct('Qrange', [64, 96], 'Lx', 15, 'Ly', 6)
        struct('Qrange', [96, 128], 'Lx', 15, 'Ly', 8)
        struct('Qrange', [128, 160], 'Lx', 15, 'Ly', 10)
        struct('Qrange', [160, 192], 'Lx', 15, 'Ly', 12)
        struct('Qrange', [192, 224], 'Lx', 15, 'Ly', 14) 
        ];
    
    materials(1).DemandFlow = demandTable.Rebar;

    model.materials.M1 = materials(1);

    % Material M2: AAC Blocks
    materials(2).ID = 'M2';
    materials(2).Name = 'AAC Blocks';
    materials(2).Unit = '1000 blocks';

    materials(2).PurchaseCost = [
        struct('Stage', [1, 2, 3], 'Quantity', [0, 10], 'PCR', 1100)
        struct('Stage', [1, 2, 3], 'Quantity', [10, 30], 'PCR', 950)
        ];

    materials(2).DeliveryCost = [
        struct('Quantity', [0, 3], 'DLC', 600)
        struct('Quantity', [3, 6], 'DLC', 1200)
        struct('Quantity', [6, 9], 'DLC', 1800)
        struct('Quantity', [9, 12], 'DLC', 2400)
        struct('Quantity', [12, 15], 'DLC', 3000)
        struct('Quantity', [15, 18], 'DLC', 3600)
        struct('Quantity', [18, 21], 'DLC', 4200)
        struct('Quantity', [21, 24], 'DLC', 4800)
        struct('Quantity', [24, 27], 'DLC', 5400)
        struct('Quantity', [27, 30], 'DLC', 6000)
        ];

    materials(2).StorageFootprint = [
        struct('Qrange', [0, 1], 'Lx', 2.5, 'Ly', 2.5)
        struct('Qrange', [1, 2], 'Lx', 5, 'Ly', 2.5)
        struct('Qrange', [2, 4], 'Lx', 5, 'Ly', 5)
        struct('Qrange', [4, 6], 'Lx', 7.5, 'Ly', 5)
        struct('Qrange', [6, 9], 'Lx', 7.5, 'Ly', 7.5)
        struct('Qrange', [9, 12], 'Lx', 10, 'Ly', 7.5)
        struct('Qrange', [12, 16], 'Lx', 10, 'Ly', 10)
        struct('Qrange', [16, 20], 'Lx', 12.5, 'Ly', 10)
        struct('Qrange', [20, 25], 'Lx', 12.5, 'Ly', 12.5)
        struct('Qrange', [25, 30], 'Lx', 15, 'Ly', 12.5)
        ];

    materials(2).DemandFlow = demandTable.AAC;
    
    model.materials.M2 = materials(2);

    % Material M3: Curtain Wall
    materials(3).ID = 'M3';
    materials(3).Name = 'Curtain Wall';
    materials(3).Unit = 'm2';

    materials(3).PurchaseCost = [
        struct('Stage', [1, 2, 3], 'Quantity', [0, 1500], 'PCR', 210)
        ];

    materials(3).DeliveryCost = [
        struct('Quantity', [0, 250], 'DLC', 1000)
        struct('Quantity', [250, 500], 'DLC', 2000)
        struct('Quantity', [500, 750], 'DLC', 3000)
        struct('Quantity', [750, 1000], 'DLC', 4000)
        struct('Quantity', [1000, 1250], 'DLC', 5000)
        struct('Quantity', [1250, 1500], 'DLC', 6000)
        ];


    materials(3).StorageFootprint = [
        struct('Qrange', [0, 250], 'Lx', 5, 'Ly', 5)
        struct('Qrange', [250, 500], 'Lx', 10, 'Ly', 5)
        struct('Qrange', [500, 750], 'Lx', 10, 'Ly', 10)
        struct('Qrange', [750, 1000], 'Lx', 15, 'Ly', 10)
        struct('Qrange', [1000, 1250], 'Lx', 20, 'Ly', 10)
        struct('Qrange', [1250, 1500], 'Lx', 20, 'Ly', 15)
        ];

    materials(3).DemandFlow = demandTable.Curtain;

    model.materials.M3 = materials(3);

    % Table 4 Materials On-Site Handling Quantities and Cost Data

    handlingData = [
    struct('MaterialID', 'M1', 'Name', 'Rebar', 'Unit', 'Ton', ...
           'Facility', 'B1', 'Stage', 1, 'RequiredQty', 286.5, ...
           'HandlingEquipment', 'Tower crane', 'HandlingQty', 2, ...
           'Speed', 5000, 'HourlyCost', 200, 'TravelCostRate',[])

    struct('MaterialID', 'M1', 'Name', 'Rebar', 'Unit', 'Ton', ...
           'Facility', 'B1', 'Stage', 2, 'RequiredQty', 280, ...
           'HandlingEquipment', 'Tower crane', 'HandlingQty', 2, ...
           'Speed', 5000, 'HourlyCost', 200, 'TravelCostRate',[])

    struct('MaterialID', 'M1', 'Name', 'Rebar', 'Unit', 'Ton', ...
           'Facility', 'B1', 'Stage', 3, 'RequiredQty', 88, ...
           'HandlingEquipment', 'Tower crane', 'HandlingQty', 2, ...
           'Speed', 5000, 'HourlyCost', 200, 'TravelCostRate',[])

    struct('MaterialID', 'M1', 'Name', 'Rebar', 'Unit', 'Ton', ...
           'Facility', 'B2', 'Stage', 3, 'RequiredQty', 124, ...
           'HandlingEquipment', 'Tower crane', 'HandlingQty', 2, ...
           'Speed', 5000, 'HourlyCost', 200, 'TravelCostRate',[])

    struct('MaterialID', 'M2', 'Name', 'AAC Blocks', 'Unit', 'M', ...
           'Facility', 'B1', 'Stage', 2, 'RequiredQty', 97.5, ...
           'HandlingEquipment', 'Tower crane', 'HandlingQty', 0.5, ...
           'Speed', 5000, 'HourlyCost', 200, 'TravelCostRate',[])

    struct('MaterialID', 'M2', 'Name', 'AAC Blocks', 'Unit', 'M', ...
           'Facility', 'B1', 'Stage', 3, 'RequiredQty', 30.5, ...
           'HandlingEquipment', 'Tower crane', 'HandlingQty', 0.5, ...
           'Speed', 5000, 'HourlyCost', 200, 'TravelCostRate',[])

    struct('MaterialID', 'M2', 'Name', 'AAC Blocks', 'Unit', 'M', ...
           'Facility', 'B2', 'Stage', 3, 'RequiredQty', 17.83, ...
           'HandlingEquipment', 'Tower crane', 'HandlingQty', 0.5, ...
           'Speed', 5000, 'HourlyCost', 200, 'TravelCostRate',[])

    struct('MaterialID', 'M3', 'Name', 'Curtain Wall', 'Unit', 'm2', ...
           'Facility', 'B1', 'Stage', 2, 'RequiredQty', 4541.67, ...
           'HandlingEquipment', 'Tower crane', 'HandlingQty', 4.5, ...
           'Speed', 5000, 'HourlyCost', 200, 'TravelCostRate',[])

    struct('MaterialID', 'M3', 'Name', 'Curtain Wall', 'Unit', 'm2', ...
           'Facility', 'B1', 'Stage', 3, 'RequiredQty', 1998.33, ...
           'HandlingEquipment', 'Tower crane', 'HandlingQty', 4.5, ...
           'Speed', 5000, 'HourlyCost', 200, 'TravelCostRate',[])

    struct('MaterialID', 'M3', 'Name', 'Curtain Wall', 'Unit', 'm2', ...
           'Facility', 'B2', 'Stage', 3, 'RequiredQty', 140, ...
           'HandlingEquipment', 'Tower crane', 'HandlingQty', 4.5, ...
           'Speed', 5000, 'HourlyCost', 200, 'TravelCostRate',[])
    ];

    model.handlingData=handlingData;

    % Table 5 Site Constraints

    % --- Safety constraints (Min distances) ---
    safetyMinPairs = {
        'B1','F2',5; 'B1','F3',5; 'B1','F7',5;
        'B2','F2',5; 'B2','F3',5; 'B2','F7',5;
        'F1','F2',15; 'F1','F3',15; 'F1','F7',15
    };

    % --- Operational constraints (Max distances) ---
    operationalMaxPairs = {
        'F1','F2',30; 'F1','F3',30; 'F1','B1',30; 'F1','B2',30;
        'F1','F4',30; 'F1','F6',30; 'F1','M1',30; 'F1','M2',30; 'F1','M3',30;
        };
    
    % --- Operational constraints (Min distances) ---
    operationalMinPairs = {
        'M3','M1',5; 'M3','M2',5; 'M3','B1',5; 'M3','B2',5
    };

    Constraints.safetyMinPairs=safetyMinPairs;
    Constraints.operationalMaxPairs=operationalMaxPairs;
    Constraints.operationalMinPairs=operationalMinPairs;

    exclusion.x = 0;
    exclusion.y = 0;
    exclusion.Lx = 15;
    exclusion.Ly = 15;
    exclusion.activeStages = [1, 2, 3];

    Constraints.exclusion=exclusion;

    model.Constraints=Constraints;

    %% 5S Parameters

    model.s5_params = struct();

    % --- وزن‌های اصلی برای تجمیع (رابطه 8) ---
    
    model.s5_params.Ws = 0.25; % وزن Sort
    model.s5_params.Wo = 0.25; % وزن SetInOrder
    model.s5_params.Wc = 0.25; % وزن Shine
    model.s5_params.Wt = 0.25; % وزن Standardize

    % --- وزن‌های داخلی (سطح 1) برای هر اصل ---
    model.s5_params.ws1 = 0.5; % وزن بخش اول P_Sort (جریمه انبارداری زودهنگام)
    model.s5_params.ws2 = 0.5; % وزن بخش دوم P_Sort (جریمه فضای اضافه)

    model.s5_params.wo1 = 0.7; % وزن بخش اول (همجواری تسهیلات)
    model.s5_params.wo2 = 0.3; % وزن بخش دوم (دسترسی به مصالح)

    model.s5_params.wc1 = 0.2; % وزن SCI
    model.s5_params.wc2 = 0.4; % وزن تداخل مسیر
    model.s5_params.wc3 = 0.4; % وزن نقض ایمنی

    % --- ضریب مقیاس آلفا برای تابع هدف نهایی ---
    % این عدد یک تصمیم مدیریتی است (مثلا 20% هزینه میانگین لجستیک)
    % فعلا یک مقدار اولیه قرار می‌دهیمSS
    model.s5_params.alpha = 1e5; 

    % --- پارامترهای SetInOrder ---
    % ماتریس وابستگی (Closeness Matrix) بر اساس تسهیلات F1 تا F7
    % ترتیب ستون و سطر: F1, F2, F3, F4, F5, F6, F7
    % مقادیر بزرگتر یعنی نزدیکی دو تسهیلات نامطلوب‌تر است
    R_matrix = [
    %   F1 F2 F3 F4 F5 F6 F7 F8
        0, 5, 5, 3, 3, 1, 5, 5;  % F1 (Crane)
        5, 0, 1, 3, 1, 3, 3, 5;  % F2 (ُSite-Offic Trailor)
        5, 1, 0, 3, 1, 3, 3, 3;  % F3 (Site-Office Trailor)
        3, 3, 3, 0, 5, 5, 1, 5;  % F4 (Fabrication Area)
        3, 1, 1, 5, 0, 3, 1, 3;  % F5 (Dump Area)
        1, 3, 3, 5, 3, 0, 1, 5;  % F6 (Lay-Down Area)
        5, 3, 3, 1, 1, 1, 0, 5;  % F7 (Labor rest Area)
        5, 5, 3, 5, 3, 5, 5, 0;  % F8 (Fuel Storage)
        ];
    model.s5_params.ClosenessMatrix = R_matrix;

    % --- پارامترهای Shine (برای محاسبه SCI) ---
    % محاسبه فضای قابل استفاده یک بار برای همیشه

    Area_site = model.site_width * model.site_height;
    Area_fixed = 0; % مقدار اولیه برای مساحت کل
    facilities_names = fieldnames(model.Facilities); % گرفتن نام تمام فسیلیتی‌ها
    
    for i = 1:length(facilities_names)
        facility_name = facilities_names{i}; % نام فسیلیتی فعلی در حلقه
        
        % بررسی اینکه آیا نام با 'B' شروع می‌شود یا نه
        if startsWith(facility_name, 'B')
            % دسترسی به فسیلیتی فعلی و محاسبه مساحت آن
            current_facility = model.Facilities.(facility_name);
            Area_fixed = Area_fixed + current_facility.Lx * current_facility.Ly;
        end
    end
    model.s5_params.UsableArea = Area_site - Area_fixed;

    % در فایل CreateModel.m، داخل model.s5_params

    % --- پارامترهای P_Shine (بخش 2 و 3) ---
    
    % تعریف کریدورهای تردد اصلی (Main Traffic Corridors)
    % هر کریدور یک مستطیل است: [x, y, Lx, Ly]
    model.s5_params.TrafficCorridors = [
        5, 5, 5, 65;   % یک کریدور عمودی در سمت چپ سایت
        5, 5, 50, 5    % یک کریدور افقی در پایین سایت
    ];
    
    % تعریف مناطق ایمنی (Safety Zones)
    % هر منطقه یک ساختار با نوع (type) و پارامترهای مشخص است
    % % مثال: محدوده چرخش جرثقیل (F1)
    % model.s5_params.SafetyZones(1).type = 'circle';
    % model.s5_params.SafetyZones(1).center_facility = 'F1'; % مرکز آن، مرکز F1 است
    % model.s5_params.SafetyZones(1).radius = 5; % شعاع چرخش 15 متر
    
    % مثال: منطقه ممنوعه نزدیک مواد قابل اشتعال (F8)
    model.s5_params.SafetyZones(1).type = 'rectangle';
    model.s5_params.SafetyZones(1).center_facility = 'F8';
    model.s5_params.SafetyZones(1).buffer_x = 2; % 5 متر حائل از هر طرف
    model.s5_params.SafetyZones(1).buffer_y = 2;

end

