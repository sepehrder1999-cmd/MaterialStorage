function [child1, child2, feasible1, feasible2] = Crossover(Parent1, Parent2, model)

    T = numel(model.Stage);
    numMaterials = length(fieldnames(model.materials));
    
    % --------- 1. تعداد FOP genes (first block) ---------
    numFOPGenes = 0;
    for t = 1:T
        for m = 1:numMaterials
            mname = ['M', num2str(m)];
            if isfield(model.materials, mname)
                demand = model.materials.(mname).DemandFlow;
                stageStart = sum([model.Stage(1:t-1).Duration]) + 1;
                stageEnd = stageStart + model.Stage(t).Duration;
                if any(demand(stageStart:min(stageEnd, end)))
                    numFOPGenes = numFOPGenes + 1;
                end
            end
        end
    end
    
    % --------- 2. Crossover: FOP Block (1-point) ---------
    totalLength = length(Parent1);
    child1 = zeros(1, totalLength);
    child2 = zeros(1, totalLength);

    cutPoint = randi([1, numFOPGenes - 1]);
    child1(1:cutPoint) = Parent1(1:cutPoint);
    child1(cutPoint+1:numFOPGenes) = Parent2(cutPoint+1:numFOPGenes);
    
    child2(1:cutPoint) = Parent2(1:cutPoint);
    child2(cutPoint+1:numFOPGenes) = Parent1(cutPoint+1:numFOPGenes);  

    % ==== 2. آپدیت Lx, Ly متریال‌ها برای هر بچه بر اساس FOP ====
    matNames = fieldnames(model.materials);

    layoutStart = numFOPGenes + 1;  % ابتدای layout در کروموزوم
    offset = 0;  % فاصله از ابتدای layout بر

    nMatTotal = 0;
    
    for t = 1:T
        stageStart = sum([model.Stage(1:t-1).Duration]) + 1;
        stageEnd = stageStart + model.Stage(t).Duration;

         % ابتدا تعداد مواد فعال در این استیج را بشماریم
        nMat = 0;
        matMap = [];  % نگه‌داری mapping شماره ماده‌ها به ایندکس داخل stage
    
        for m = 1:length(matNames)
            mname = matNames{m};
            demandFlow = model.materials.(mname).DemandFlow;
    
            if length(demandFlow) < stageEnd
                relevantDemand = demandFlow(stageStart:end);
            else
                relevantDemand = demandFlow(stageStart:stageEnd);
            end
    
            if isempty(relevantDemand) || all(relevantDemand == 0)
                continue;
            end
    
            % === مقدار FOP جدید را بخوان از ژن های کودک ===
            nMatTotal = nMatTotal + 1;
            nMat = nMat + 1;
            matMap(end+1) = m;  % نگه‌داری mapping

            fop1 = child1(nMatTotal);
            fop2 = child2(nMatTotal);
    
            % ===== محاسبه Lx, Ly بر اساس FOP برای هر کودک =====
            [lx1, ly1] = GetFootprintFromFOP(model.materials.(mname), fop1, relevantDemand, t, model);
            [lx2, ly2] = GetFootprintFromFOP(model.materials.(mname), fop2, relevantDemand, t, model);
    
            % موقعیت دقیق در کروموزوم:
            idx1 = layoutStart + offset + (nMat - 1) * 5;

            % ابتدا مقدار theta فعلی هر بچه را بخوان
            theta1 = Parent1(idx1 + 2);
            theta2 = Parent2(idx1 + 2);
            
            % اگر 90 درجه بودن، جای lx و ly را عوض کن
            if mod(theta1, 180) == 90
                [lx1, ly1] = deal(ly1, lx1);
            end
            if mod(theta2, 180) == 90
                [lx2, ly2] = deal(ly2, lx2);
            end

            child1(idx1 + 3) = lx1;
            child1(idx1 + 4) = ly1;
    
            child2(idx1 + 3) = lx2;
            child2(idx1 + 4) = ly2;
        end

        % حالا تعداد facilityهای این stage را حساب کن تا offset بعدی مشخص شود
        nFac = 0;
        allFacNames  = fieldnames(model.Facilities);
        facNames = allFacNames(startsWith(allFacNames, 'F'));
        for f = 1:length(facNames)
            fObj = model.Facilities.(facNames{f});
            if ismember(t, fObj.activestages) && startsWith(facNames{f}, 'F')
                nFac = nFac + 1;
            end
        end

        % برای مرحله بعدی: offset را افزایش بده (برای کل stage فعلی)
        offset = offset + 5 * (nMat + nFac);
    end

    % --- 3. Layout Block Crossover ---
    pos1 = numFOPGenes + 1;
    
    for t = 1:T
        nMat = 0;
        stageStart = sum([model.Stage(1:t-1).Duration]) + 1;
        stageEnd = stageStart + model.Stage(t).Duration;
    
        % --- شمارش مواد ---
        for m = 1:length(matNames)
            mname = matNames{m};
            if any(model.materials.(mname).DemandFlow(stageStart:min(stageEnd, end)))
                nMat = nMat + 1;
            end
        end
    
        % --- استخراج Facilityهای این stage ---
        activeFacs = {};
        stationaryIdx = [];
        moveableIdx = [];
    
        for f = 1:length(facNames)
            fname = facNames{f};
            fobj = model.Facilities.(fname);
            
            % فقط اگر در این استیج فعال است:
            if ismember(t, fobj.activestages)
                activeFacs{end+1} = fname;
                
                % فقط فسیلیتی‌هایی که با 'F' شروع می‌شوند بررسی شوند
                if isfield(fobj, 'Type') && strcmpi(fobj.Type, 'S')
                    stationaryIdx(end+1) = length(activeFacs);  % index نسبی در activeFacs
                else
                    moveableIdx(end+1) = length(activeFacs);
                end
            end
        end
    
        % کل آیتم‌ها
        numItems = nMat + length(activeFacs);
        blockSize = 5 * numItems;
    
        % تعیین نوع هر آیتم در stage فعلی
        itemTypes = cell(numItems, 1);  % 'M' یا 'F'
        itemNames = cell(numItems, 1); % برای شناسایی دقیق‌تر
        idxItem = 1;
    
        % ابتدا مواد
        for m = 1:length(matNames)
            mname = matNames{m};
            if any(model.materials.(mname).DemandFlow(stageStart:min(stageEnd, end)))
                itemTypes{idxItem} = 'M';
                itemNames{idxItem} = mname;
                idxItem = idxItem + 1;
            end
        end
    
        % سپس فسیلیتی‌ها
        for i = 1:length(activeFacs)
            itemTypes{idxItem} = 'F';
            itemNames{idxItem} = activeFacs{i};
            idxItem = idxItem + 1;
        end
    
        % --- اعمال Crossover برای هر آیتم ---
        cutPoint = randi([1, numItems - 1]);
    
        for i = 1:numItems
            idx = pos1 + (i - 1) * 5;
                
            if strcmp(itemTypes{i}, 'F')
                facName = itemNames{i};
                
                % فقط اگر اسم با 'F' شروع شود بررسی شود
                if startsWith(facName, 'F')
                    fobj = model.Facilities.(facName);
                    
                    % بررسی ایمن وجود فیلد Type
                    if isfield(fobj, 'Type') && strcmpi(fobj.Type, 'S')
                        % کپی مستقیم بدون تغییر
                        child1(idx:idx+4) = Parent1(idx:idx+4);
                        child2(idx:idx+4) = Parent2(idx:idx+4);
                    else
                        % مواد یا Facilities متحرک: crossover معمولی
                        if i <= cutPoint
                            child1(idx:idx+2) = Parent1(idx:idx+2);
                            child2(idx:idx+2) = Parent2(idx:idx+2);
                        else
                            child1(idx:idx+2) = Parent2(idx:idx+2);
                            child2(idx:idx+2) = Parent1(idx:idx+2);
                        end
                    
                        % --- Lx, Ly ---
                        if strcmp(itemTypes{i}, 'M')  % برای مواد
                            child1(idx+3:idx+4) = child1(idx+3:idx+4);
                            child2(idx+3:idx+4) = child2(idx+3:idx+4);
                        else  % برای Facilities متحرک
                            if i <= cutPoint
                                child1(idx+3:idx+4) = Parent1(idx+3:idx+4);
                                child2(idx+3:idx+4) = Parent2(idx+3:idx+4);
                            else
                                child1(idx+3:idx+4) = Parent2(idx+3:idx+4);
                                child2(idx+3:idx+4) = Parent1(idx+3:idx+4);
                            end
                        end
                    end
                end
            end
    
        %     % مواد یا فسیلیتی‌های moveable: crossover معمولی
        %     if i <= cutPoint
        %         child1(idx:idx+2) = Parent1(idx:idx+2);
        %         child2(idx:idx+2) = Parent2(idx:idx+2);
        %     else
        %         child1(idx:idx+2) = Parent2(idx:idx+2);
        %         child2(idx:idx+2) = Parent1(idx:idx+2);
        %     end
        % 
        %     % --- Lx, Ly ---
        %     if strcmp(itemTypes{i}, 'M')  % برای مواد، نگه دار
        %         % child1(idx+3:idx+4) قبلاً set شده
        %         child1(idx+3:idx+4) = child1(idx+3:idx+4);
        %         child2(idx+3:idx+4) = child2(idx+3:idx+4);
        %     else  % برای فسیلیتی‌های moveable، crossover
        %         if i <= cutPoint
        %             child1(idx+3:idx+4) = Parent1(idx+3:idx+4);
        %             child2(idx+3:idx+4) = Parent2(idx+3:idx+4);
        %         else
        %             child1(idx+3:idx+4) = Parent2(idx+3:idx+4);
        %             child2(idx+3:idx+4) = Parent1(idx+3:idx+4);
        %         end
        %     end
        end
    
        pos1 = pos1 + blockSize;
    end

    [child1, feasible1] = RepairChild(child1, model);
    [child2, feasible2] = RepairChild(child2, model);

end



function [Lx, Ly] = GetFootprintFromFOP(material, fop, demand, stageIndex, model)

demandFlow = material.DemandFlow;
            stages = model.Stage;
            stageDuration = stages(stageIndex).Duration;
    if stageIndex == 1
       stageStart = 1;
       stageEnd = stageStart + stageDuration;

    else
       stageStart = sum([stages(1:stageIndex-1).Duration]) + 1;
       stageEnd = stageStart + stageDuration;
    end

                    
    if length(demandFlow) < stageEnd
       relevantDemand = demandFlow(stageStart:end);
    else
       relevantDemand = demandFlow(stageStart:stageEnd);
    end

    if length(demand) < fop
        orderQty = demand(end);
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
                w2 = s * fop + 1;
    
                if w1 > length(demand)
                    continue;
                end
                if w2 > length(demand)
                    w2 = length(demand);
                end
    
                if stageIndex == 1 && s == 1
                    orderQty(s) = demand(w2);
                    continue
                end
    
                orderQty(s) = demand(w2) - demand(w1);
            end
        end
    end

    maxQty = max(orderQty);
    % sol.Stage(i).Materials.(mat).orderQty = orderQty;

    
    footprintOptions = material.StorageFootprint;

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
end

function [child_fixed, feasible] = RepairChild(child, model)

    sol = ParseSolution(child, model);  % مرحله 1: پارس کردن

    T = numel(model.Stage);
    gridW = model.site_width;
    gridH = model.site_height;
    feasible = true;

    % === مرحله به مرحله
    for t = 1:T
        gridMap = zeros(gridH, gridW);  % مرحله 2: ساخت گرید

        % مرحله 3: تعریف blockedZones
        blockedZones = struct('x', {}, 'y', {}, 'width', {}, 'height', {});
        ex = model.Constraints.exclusion;
        blockedZones(end+1) = struct('x', ex.x, 'y', ex.y, 'width', ex.Lx, 'height', ex.Ly);

        % مرحله 4: قرار دادن B و G ها روی سایت
        if isfield(sol.Stage(t), 'Buildings')
            all_bnames = fieldnames(sol.Stage(t).Buildings);
            bnames = all_bnames(startsWith(all_bnames, 'B'));
            for i = 1:length(bnames)
                b = sol.Stage(t).Buildings.(bnames{i});
                bx = max(1, floor(b.x));
                by = max(1, floor(b.y));
                bw = floor(b.Lx);
                bh = floor(b.Ly);
                gridMap(by:by+bh-1, bx:bx+bw-1) = 1;

                blockedZones(end+1) = struct('x', bx, 'y', by, 'width', bw, 'height', bh);
            end
        end

        if isfield(sol.Stage(t), 'Facilities')
            facNames = fieldnames(sol.Stage(t).Facilities);
            for i = 1:length(facNames)
                fname = facNames{i};
                if isfield(model.Facilities.(fname), 'Type') && strcmpi(model.Facilities.(fname).Type, 'S')
                    f = sol.Stage(t).Facilities.(fname);
                    fx = max(1, floor(f.x));
                    fy = max(1, floor(f.y));
                    fw = floor(f.Lx);
                    fh = floor(f.Ly);
        
                    blockedZones(end+1) = struct('x', fx, 'y', fy, 'width', fw, 'height', fh);
                end
            end
        end


        % === مرتب‌سازی فسیلیتی‌ها بر اساس مساحت نزولی ===
        if isfield(sol.Stage(t), 'Facilities')
            facNames = fieldnames(sol.Stage(t).Facilities);
        
            % فقط فسیلیتی‌هایی که در مدل تعریف شده‌اند (برای اطمینان)
            areas = zeros(size(facNames));
            for j = 1:length(facNames)
                fname = facNames{j};
                if isfield(model.Facilities, fname)
                    f = model.Facilities.(fname);
                    if isfield(f, 'Lx') && isfield(f, 'Ly')
                        areas(j) = f.Lx * f.Ly;
                    else
                        areas(j) = 0; % اگر مشخص نشده بود
                    end
                end
            end
        
            [~, sortIdx] = sort(areas, 'descend');
            facNames = facNames(sortIdx);  % facNames حالا مرتب شده است
        end




        % مرحله 5: قرار دادن فسیلیتی‌های ثابت
        if isfield(sol.Stage(t), 'Facilities')
            for i = 1:length(facNames)
                fname = facNames{i};
                f = sol.Stage(t).Facilities.(fname);
                if isfield(model.Facilities.(fname), 'Type') && strcmpi(model.Facilities.(fname).Type, 'S')
                    % if ~checkFree(gridMap, f.x, f.y, f.Lx, f.Ly, blockedZones)
                    %     feasible = false;
                    %     break;
                    % end
                    gridMap(f.y:f.y+f.Ly-1, f.x:f.x+f.Lx-1) = 1;
                end
            end
        end

        % مرحله 6: قرار دادن فسیلیتی‌های متحرک با randomWalk
        for i = 1:length(facNames)
            fname = facNames{i};
            f = sol.Stage(t).Facilities.(fname);
            if isfield(model.Facilities.(fname), 'Type') && strcmpi(model.Facilities.(fname).Type, 'M')
                for tries = 1:10000
                    if checkFree(gridMap, f.x, f.y, f.Lx, f.Ly, blockedZones)
                        gridMap(f.y:f.y+f.Ly-1, f.x:f.x+f.Lx-1) = 1;
                        sol.Stage(t).Facilities.(fname) = f;
                        break;
                    else
                        [f.x, f.y] = randomWalk(f.x, f.y, gridW, gridH, f.Lx, f.Ly);
                    end
                end
                if tries == 10000
                    feasible = false;
                    break;
                end
            end
        end

        % مرحله 7: جانمایی متریال‌ها
        if isfield(sol.Stage(t), 'Materials')
            matNames = fieldnames(sol.Stage(t).Materials);
            Placed = false;

            while ~Placed
                 allPlaced = true;

                for i = 1:length(matNames)
                    mname = matNames{i};
                    m = sol.Stage(t).Materials.(mname);
                     placed_this_material = false;
    
                    for tries = 1:10000
                        if checkFree(gridMap, m.x, m.y, m.Lx, m.Ly, blockedZones)
                            gridMap(m.y:m.y+m.Ly-1, m.x:m.x+m.Lx-1) = 1;
                            sol.Stage(t).Materials.(mname) = m;
                             placed_this_material = true;
                            break;
                        else
                            [m.x, m.y] = randomWalk(m.x, m.y, gridW, gridH, m.Lx, m.Ly);
                        end
                    end
                    if ~placed_this_material
                        mat = model.materials.(mname);
                        stageStart = sum([model.Stage(1:t-1).Duration]) + 1;
                        stageEnd = stageStart + model.Stage(t).Duration;
                        relevantDemand = mat.DemandFlow(stageStart:min(stageEnd, end));
                        fop0 = m.FOP;
                        fopOptions = (fop0 - 1):-1:1;

                        for fop = fopOptions
                            % محاسبه maxQty
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
                                        w2 = s * fop + 1;
                
                                        if w1 > length(relevantDemand)
                                            continue;
                                        end
                
                                        if w2 > length(relevantDemand)
                                            w2 = length(relevantDemand);
                                        end
                
                                        if i == 1 && s == 1
                                            orderQty(s) = relevantDemand(w2);
                                            continue;
                                        end
                
                                        orderQty(s) = relevantDemand(w2) - relevantDemand(w1);
                                    end
                                end
                            end
                            maxQty = max(orderQty);
                            sol.Stage(i).Materials.(mname).orderQty = orderQty;

    
                            % ابعاد جدید بر اساس StorageFootprint
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
            
                            x = floor(m.x);
                            y = floor(m.y);

                            placed_this_material = false;
    
                            for tries = 1:10000
                                if checkFree(gridMap, x, y, Lx, Ly, blockedZones)
                                    gridMap(y:y+Ly-1, x:x+Lx-1) = 1;
                                    m.x=x;
                                    m.y=y;
                                    m.Lx=Lx;
                                    m.Ly=Ly;
                                    m.FOP=fop;
                                    sol.Stage(t).Materials.(mname) = m;
                                    
                                    placed_this_material = true;
                                    break;
                                else
                                    [x, y] = randomWalk(x, y, gridW, gridH, Lx, Ly);
                                end
                            end
                            if placed_this_material
                                 break;
                            % elseif ~placed_this_material
                            end
                        end
                    end
            % اگر بعد از کاهش FOP هم جا نشد → infeasible
                    if ~placed_this_material
                        allPlaced = false;
                        break;
                    end
                end
                % اگر همه متریال‌ها موفق شدن، while تموم شه
                if allPlaced
                    Placed = true;
                else
                    feasible = false;
                    break;
                end
            end
        end
    end

    % مرحله 8: بازسازی ژنوم از روی sol
    sol.pos = [];
    % matNames = fieldnames(model.materials);

    % --- FOP ها
    for t = 1:T
        matNames = fieldnames(sol.Stage(t).Materials);
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

    % --- Layout: مواد و فسیلیتی‌ها
    for t = 1:T
        % مواد
        for m = 1:length(matNames)
            mname = matNames{m};
            if isfield(sol.Stage(t).Materials, mname)
                mat = sol.Stage(t).Materials.(mname);
                if all(isfield(mat, {'x', 'y', 'theta', 'Lx', 'Ly'}))
                    sol.pos(end+1:end+5) = [mat.x, mat.y, mat.theta, mat.Lx, mat.Ly];
                end
            end
        end

        % فسیلیتی‌ها (مرتب شده)
        if isfield(sol.Stage(t), 'Facilities')
            facNames = fieldnames(sol.Stage(t).Facilities);
            fFiltered = facNames(startsWith(facNames, 'F'));
            fNumbers = cellfun(@(s) sscanf(s, 'F%d'), fFiltered);
            [~, sortIdx] = sort(fNumbers);
            sortedF = fFiltered(sortIdx);
            for j = 1:length(sortedF)
                fname = sortedF{j};
                f = sol.Stage(t).Facilities.(fname);
                if all(isfield(f, {'x', 'y', 'theta', 'Lx', 'Ly'}))
                    sol.pos(end+1:end+5) = [f.x, f.y, f.theta, f.Lx, f.Ly];
                end
            end
        end
    end

    child_fixed = sol.pos;

end
