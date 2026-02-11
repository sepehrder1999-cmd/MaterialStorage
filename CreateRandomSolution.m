function sol = CreateRandomSolution(model)

    T = numel(model.Stage);
    [H, W] = deal(model.site_height, model.site_width);
    
    blockedZones = struct('x', {}, 'y', {}, 'width', {}, 'height', {});
    
    % --- Exclusion Zone ---
    blockedZones(end+1) = struct( ...
        'x', model.Constraints.exclusion.x, ...
        'y', model.Constraints.exclusion.y, ...
        'width', model.Constraints.exclusion.Lx, ...
        'height', model.Constraints.exclusion.Ly);

    nBlockedBefore = numel(blockedZones);

    
    % --- Buildings (with Lx and Ly → mapped to width, height) ---
    facilityNames = fieldnames(model.Facilities);
    isBuilding = startsWith(facilityNames, 'B');   
    buildingNames = facilityNames(isBuilding);    

    for i = 1:T
            gridMap = zeros(H, W);
            sol.Stage(i).occupiedRects = struct( ...
                'name', {}, ...
                'x', {}, ...
                'y', {}, ...
                'width', {}, ...
                'height', {});
       
                    % لیست فسیلیتی فعال
            activeFacs = {};
            for j = 1:length(facilityNames)
                fname = facilityNames{j};
                f = model.Facilities.(fname);
                if ismember(i, f.activestages)
                    activeFacs{end+1} = fname;
                end
            end
        
            
            % ساختمان‌های ثابت (B1، B2)
            fname = startsWith(facilityNames, {'B'});
            fname = facilityNames(fname);  % فقط B و Gها
            for j = 1:numel(fname)
                name = fname{j};  % مقدار رشته‌ای از cell array
                
                f = model.Facilities.(name);
                
                if ismember(i, f.activestages)
                    x = f.x; y = f.y; w = f.Lx; h = f.Ly;
                    if x == 0 || y == 0
                        gridMap(y+1:y+h, x+1:x+w) = 1;
                    else
                        gridMap(y:y+h-1, x:x+w-1) = 1;
                    end
            
                    sol.Stage(i).Buildings.(name) = struct('x',x,'y',y,'Lx',w,'Ly',h);
                    sol.Stage(i).occupiedRects(end+1) = struct( ...
                               'name', name, ...
                               'x', x, ...
                               'y', y, ...
                               'width', w, ...
                               'height', h);
                end
            end


            
            % مرتب‌سازی بر اساس مساحت (نزولی)
            areas = zeros(size(activeFacs));
            for j = 1:length(activeFacs)
                f = model.Facilities.(activeFacs{j});
                areas(j) = f.Lx * f.Ly;
            end
            [~, sortIdx] = sort(areas, 'descend');
            sortedFacs = activeFacs(sortIdx);
        
            thetaOptions = [0, 90];
            % چرخش ابعاد در صورت نیاز
            for j = 1:length(sortedFacs)
                fname = sortedFacs{j};
        
                if startsWith(fname, 'F')
                    theta = thetaOptions(randi([1, 2]));
                    model.Facilities.(fname).theta = theta;
                    if theta == 90
                        % جابجایی ابعاد
                        temp = model.Facilities.(fname).Lx;
                        model.Facilities.(fname).Lx = model.Facilities.(fname).Ly;
                        model.Facilities.(fname).Ly = temp;
                        model.Facilities.(fname).theta = theta;
                    end
                end
            end
            
            % === مرحله اول: جانمایی فسیلیتی‌ ها ===
            for j = 1:length(sortedFacs)
                fname = sortedFacs{j};
                if startsWith(fname, 'F')
        
                    f = model.Facilities.(fname);
            
                    if isfield(f, 'Type') && strcmp(f.Type, 'S')
                        % اگر قبلاً در مرحله 1 قرار گرفته، استفاده کن
                        if i == 1
                            placed = false; trial = 0;
                            x = randi([1, W - f.Lx + 1]);
                            y = randi([1, H - f.Ly + 1]);

                            for c = 1:numel(buildingNames)
                                b = model.Facilities.(buildingNames{c});
                                if all(isfield(b, {'x', 'y', 'Lx', 'Ly'}))
                                    blockedZones(end+1) = struct( ...
                                        'x', b.x, ...
                                        'y', b.y, ...
                                        'width', b.Lx, ...
                                        'height', b.Ly);
                                end
                            end

                            while ~placed && trial < 10000
                                trial = trial + 1;

                                if checkFree(gridMap, x, y, f.Lx, f.Ly, blockedZones)
                                    gridMap(y:y+f.Ly-1, x:x+f.Lx-1) = 1;
                                    sol.Stage(i).Facilities.(fname) = struct('x',x,'y',y,'Lx',f.Lx,'Ly',f.Ly,'theta', model.Facilities.(fname).theta);
                                    sol.Stage(i).occupiedRects(end+1) = struct( ...
                                                'name', fname, ...
                                                'x', x, ...
                                                'y', y, ...
                                                'width', f.Lx, ...
                                                'height', f.Ly);
                                    placed = true;
                                else
                                    [x, y] = randomWalk(x, y, W, H, f.Lx, f.Ly);
                                    % [x, y] = FindFreePositionGreedy(gridMap, f.Lx, f.Ly, blockedZones);
                                end
                            end
                        else
                            % از محل مرحله قبل کپی شود
                            prev = sol.Stage(i-1).Facilities.(fname);
                            sol.Stage(i).Facilities.(fname) = prev;
                            gridMap(prev.y:prev.y+prev.Ly-1, prev.x:prev.x+prev.Lx-1) = 1;
                            sol.Stage(i).occupiedRects(end+1) = struct( ...
                                    'name', fname, ...
                                    'x', prev.x, ...
                                    'y', prev.y, ...
                                    'width', prev.Lx, ...
                                    'height', prev.Ly);
                        end
                    else
                    % === مرحله دوم: بقیه فسیلیتی‌ها (قابل جابجایی) ===

                    blockedZones(nBlockedBefore+1:end) = [];
                    if isfield(sol.Stage(i), 'Buildings')
                        all_bnames = fieldnames(sol.Stage(i).Buildings);
                        bnames = all_bnames(startsWith(all_bnames, 'B'));
                        for a = 1:length(bnames)
                            b = sol.Stage(i).Buildings.(bnames{a});
                            bx = max(1, floor(b.x));
                            by = max(1, floor(b.y));
                            bw = floor(b.Lx);
                            bh = floor(b.Ly);
    
                            blockedZones(end+1) = struct('x', bx, 'y', by, 'width', bw, 'height', bh);
                        end
                        if i>1
                           allfnames = fieldnames(model.Facilities);
                            fnames = allfnames(startsWith(allfnames, 'F'));
                            for a = 1:length(fnames)
                                ff = model.Facilities.(fnames{a});
                                if isfield(ff, 'Type') && strcmp(ff.Type, 'S')
                                    fx = sol.Stage(i-1).Facilities.(fnames{a}).x;
                                    fy = sol.Stage(i-1).Facilities.(fnames{a}).y;
                                    fw = sol.Stage(i-1).Facilities.(fnames{a}).Lx;
                                    fh = sol.Stage(i-1).Facilities.(fnames{a}).Ly;
                            
                                    blockedZones(end+1) = struct('x', fx, 'y', fy, 'width', fw, 'height', fh);
                                end
                            end
                        end
                    end

            
                        placed = false; trial = 0;
                        x = randi([1, W - f.Lx + 1]);
                        y = randi([1, H - f.Ly + 1]);
                        while ~placed && trial < 10000
                            trial = trial + 1;
                            if checkFree(gridMap, x, y, f.Lx, f.Ly, blockedZones)
                                gridMap(y:y+f.Ly-1, x:x+f.Lx-1) = 1;
                                sol.Stage(i).Facilities.(fname) = struct('x',x,'y',y,'Lx',f.Lx,'Ly',f.Ly,'theta', model.Facilities.(fname).theta);
                                sol.Stage(i).occupiedRects(end+1) = struct( ...
                                            'name', fname, ...
                                            'x', x, ...
                                            'y', y, ...
                                            'width', f.Lx, ...
                                            'height', f.Ly);
                                placed = true;
                            else
                                [x, y] = randomWalk(x, y, W, H, f.Lx, f.Ly);
                                % [x, y] = FindFreePositionGreedy(gridMap, f.Lx, f.Ly, blockedZones);
                            end
                        end
                        if ~placed
                            % تلاش مجدد با چرخاندن facility (تغییر تتا)
                            f = model.Facilities.(fname);  % فسیلیتی را بگیر
                            originalLx = f.Lx;
                            originalLy = f.Ly;
                           
                            % تغییر تتا و ابعاد
            
                            if ~isfield(f, 'theta')
                                f.theta = 0;
                            end
                            theta = 90 - mod(f.theta, 180);
                                  % اگر 0 بود بشه 90 و برعکس
                            model.Facilities.(fname).theta = theta;
                            model.Facilities.(fname).Lx = originalLy;
                            model.Facilities.(fname).Ly = originalLx;
                            
                            % تلاش مجدد برای جایابی
                            placed = false;
                            trial = 0;
                            x = randi([1, W - model.Facilities.(fname).Lx + 1]);
                            y = randi([1, H - model.Facilities.(fname).Ly + 1]);
                        
                            while ~placed && trial < 10000
                                trial = trial + 1;
                                if checkFree(gridMap, x, y, model.Facilities.(fname).Lx, model.Facilities.(fname).Ly, blockedZones)
                                    gridMap(y:y+f.Ly-1, x:x+f.Lx-1) = 1;
                                    sol.Stage(i).Facilities.(fname) = struct('x',x,'y',y,'Lx',f.Lx,'Ly',f.Ly,'theta', model.Facilities.(fname).theta);
                                    sol.Stage(i).occupiedRects(end+1) = struct( ...
                                                'name', fname, ...
                                                'x', x, ...
                                                'y', y, ...
                                                'width', f.Lx, ...
                                                'height', f.Ly);
                                    placed = true;
                                 else
                                    [x, y] = randomWalk(x, y, W, H, model.Facilities.(fname).Lx, model.Facilities.(fname).Ly);
                                end
                            end
                            % warning('فسیلیتی %s در مرحله %d جا نشد → تلاش دوباره', fname, i);
                            sol = CreateRandomSolution(model);
                            return;
                        end
                    end
                end
            end

            % استخراج و مرتب‌سازی ساختار Facilities برای مرحله i
            facNames = fieldnames(sol.Stage(i).Facilities);
                
            % فیلتر فقط فسیلیتی‌هایی که با 'F' شروع می‌شن
            fFiltered = facNames(startsWith(facNames, 'F'));
                
            % استخراج شماره‌ی عددی فسیلیتی‌ها برای مرتب‌سازی
            fNumbers = cellfun(@(s) sscanf(s, 'F%d'), fFiltered);
                
            % مرتب‌سازی بر اساس شماره
            [~, sortIdx] = sort(fNumbers);
                
            % ساخت ساختار مرتب‌شده جدید
            sortedFacilities = struct();
            for k = 1:length(sortIdx)
                fname = fFiltered{sortIdx(k)};
                sortedFacilities.(fname) = sol.Stage(i).Facilities.(fname);
            end
                
            % جایگزینی در ساختار اصلی
            sol.Stage(i).Facilities = sortedFacilities;


        
                % === 4. جانمایی متریال‌ها ===
            stages = model.Stage;
            stageDuration = stages(i).Duration;
            if i == 1
                stageStart = 1;
                stageEnd = stageStart + stageDuration;

            else
                stageStart = sum([stages(1:i-1).Duration]) + 1;
                stageEnd = stageStart + stageDuration;

            end

            matNames = fieldnames(model.materials);
            fopOptions = [4, 3, 2, 1];  % هفته به ترتیب ترجیح

            thetaVec = zeros(length(matNames), 1);
            for t = 1:length(matNames)
                thetaVec(t) = 90 * randi([0 1]);  % 0 یا 90 به صورت تصادفی
            end

            for m = 1:length(matNames)
                fieldname = matNames{m};
                mat = model.materials.(fieldname);
                demandFlow = mat.DemandFlow;

                    if length(demandFlow) < stageEnd
                       relevantDemand = demandFlow(stageStart:end);
                    else
                        relevantDemand = demandFlow(stageStart:stageEnd);
                    end

                    if isempty(relevantDemand) || all(relevantDemand == 0)
                        continue;
                    end

                placed = false;
                for fop = fopOptions
                    % محاسبه بیشترین مقدار سفارش
                    if length(relevantDemand) < fop
                       orderQty = relevantDemand(end);
                    else
                        if fop > 1
                            nSteps = floor(length(relevantDemand) / fop) + 1;
                        elseif fop == 1
                            nSteps = floor(length(relevantDemand) / fop);
                        else
                            nSteps = length(relevantDemand);
                        end
                        orderQty = zeros(1, nSteps);
                        if fop == 0
                            if ~isempty(relevantDemand)
                                orderQty = [relevantDemand(2)-demandFlow(stageStart); diff(relevantDemand)];
                            else
                                orderQty = [];
                            end
                        else
                            for s = 1:nSteps
                                w1 = (s-1)*fop + 1;
                                w2 = s * fop+1;
    
                                if w1 > length(relevantDemand)
                                    continue;  % این بازه دیگه در دسترس نیست
                                end
    
                               if w2 > length(relevantDemand)
                                  w2 = length(relevantDemand);
                               end
    
                               % پرش از بازه‌های نامعتبر
                                if i == 1 && s == 1 && i == 1
                                   orderQty(s) = relevantDemand(w2);
                                   continue
                                end
    
                                orderQty(s) = relevantDemand(w2) - relevantDemand(w1);
                            end
                        end
                   end
                        maxQty = max(orderQty);
                        sol.Stage(i).Materials.(fieldname).orderQty = orderQty;

                    % انتخاب ابعاد ذخیره‌سازی مناسب
                    footprintOptions = mat.StorageFootprint;
                    selectedFootprint = [];
                    for f = 1:length(footprintOptions)
                        qrange = footprintOptions(f).Qrange;
                        if maxQty >= qrange(1) && maxQty <= qrange(2)
                           selectedFootprint = footprintOptions(f);
                            break;
                       end
                    end
                    if isempty(selectedFootprint)
                        selectedFootprint = footprintOptions(end);
                    end

                    Lx = floor(selectedFootprint.Lx);
                    Ly = floor(selectedFootprint.Ly);

                    theta = thetaVec(m);
                   if theta == 90
                       [Lx, Ly] = deal(Ly, Lx);
                   end

                   posX = randi([1, W - Lx + 1]);
                   posY = randi([1, H - Ly + 1]);

                   trial = 0;

                    while ~placed && trial < 10000
                        trial = trial + 1;
                        if checkFree(gridMap, posX, posY, Lx, Ly, blockedZones)
                            % مقداردهی موقت
                            sol.Stage(i).Materials.(fieldname).x = posX;
                            sol.Stage(i).Materials.(fieldname).y = posY;
                            sol.Stage(i).Materials.(fieldname).Lx = Lx;
                            sol.Stage(i).Materials.(fieldname).Ly = Ly;
                            sol.Stage(i).Materials.(fieldname).FOP = fop;
                            sol.Stage(i).Materials.(fieldname).theta = theta;


                            sol.Stage(i).occupiedRects(end+1) = struct( ...
                                        'name', fieldname, ...
                                        'x', posX, ...
                                        'y', posY, ...
                                        'width', Lx, ...
                                        'height', Ly);
                            gridMap(posY:posY+Ly-1, posX:posX+Lx-1) = 1;
                            placed = true;

                        else
                            [posX, posY] = randomWalk(posX, posY, W, H, Lx, Ly);
                        end
                    end

                    if placed
                        break;  % وقتی جانمایی موفق بود، سراغ FOP بعدی نرو
                    end

                    [Lx, Ly] = deal(Ly, Lx);
                    if theta == 90
                        theta = 0;
                    else
                        theta= 90;
                    end
                    trial = 0;
                    while ~placed & trial < 10000
                        trial = trial + 1;
                        if checkFree(gridMap, posX, posY, Lx, Ly, blockedZones)
                            % مقداردهی موقت
                            sol.Stage(i).Materials.(fieldname).x = posX;
                            sol.Stage(i).Materials.(fieldname).y = posY;
                            sol.Stage(i).Materials.(fieldname).Lx = Lx;
                            sol.Stage(i).Materials.(fieldname).Ly = Ly;
                            sol.Stage(i).Materials.(fieldname).FOP = fop;
                            sol.Stage(i).Materials.(fieldname).theta = theta;


                            sol.Stage(i).occupiedRects(end+1) = struct( ...
                                        'name', fieldname, ...
                                        'x', posX, ...
                                        'y', posY, ...
                                        'width', Lx, ...
                                        'height', Ly);
                            gridMap(posY:posY+Ly-1, posX:posX+Lx-1) = 1;
                            placed = true;

                        else
                            [posX, posY] = randomWalk(posX, posY, W, H, Lx, Ly);
                        end
                    end
                    if placed
                        break;  % وقتی جانمایی موفق بود، سراغ FOP بعدی نرو
                    end

                end

                if ~placed
                    % warning('متریال %s در مرحله %d جانمایی نشد.', fieldname, i);
                    sol = CreateRandomSolution(model);
                end
            end
            sol.Stage(i).gridMap = gridMap;
    end

    sol.pos = [];
    matNames = fieldnames(model.materials);


        % --------- 1. FOP ها به ترتیب stage و material ---------
    for t = 1:T
        for m = 1:length(matNames)
            mname = matNames{m};
            if isfield(sol.Stage(t).Materials, mname)
                mat = sol.Stage(t).Materials.(mname);
                if isfield(mat, 'FOP') && ~isempty(mat.FOP)
                    sol.pos(end+1) = mat.FOP;
                end
            end
        end
    end
    
    % --------- 2. برای هر stage: مواد و سپس facilities (x, y, theta, Lx, Ly) ---------
    for t = 1:T
        % --- 2.1: Materials ---
        for m = 1:length(matNames)
            mname = matNames{m};
            if isfield(sol.Stage(t).Materials, mname)
                mat = sol.Stage(t).Materials.(mname);
                if all(isfield(mat, {'x', 'y', 'theta', 'Lx', 'Ly'}))
                    sol.pos(end+1:end+5) = [mat.x, mat.y, mat.theta, mat.Lx, mat.Ly];
                end
            end
        end
    
        % --- 2.2: Facilities (فقط Fهایی که مرتب شدن) ---
        if isfield(sol.Stage(t), 'Facilities')
            facNames = fieldnames(sol.Stage(t).Facilities);
        
            for j = 1:length(facNames)
                fname = facNames{j};
                f = sol.Stage(t).Facilities.(fname);
                if all(isfield(f, {'x', 'y', 'theta', 'Lx', 'Ly'}))
                    sol.pos(end+1:end+5) = [f.x, f.y, f.theta, f.Lx, f.Ly];
                end
            end
        end
    end
end
