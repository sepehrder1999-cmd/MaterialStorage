function [mutated, feasible] = Mutation(chromosome, model, mu)

    % مقدارهای ممکن برای FOP
    FOP_values = [1, 2, 3, 4];

    T = numel(model.Stage);
    matNames = fieldnames(model.materials);
    facNames = fieldnames(model.Facilities);
    numMaterials = length(matNames);
    
    mutated = chromosome;  % کپی اولیه

    % --------- 1. FOP Mutation ---------
    numFOPGenes = 0;
    for t = 1:T
        for m = 1:numMaterials
            mname = matNames{m};
            demand = model.materials.(mname).DemandFlow;

            stageStart = sum([model.Stage(1:t-1).Duration]) + 1;
            stageEnd = stageStart + model.Stage(t).Duration;

            if any(demand(stageStart:min(stageEnd, end)))
                numFOPGenes = numFOPGenes + 1;

                if rand < mu
                    oldVal = mutated(numFOPGenes);
                    newVal = randsample(setdiff(FOP_values, oldVal), 1);
                    mutated(numFOPGenes) = newVal;
                end
            end
        end
    end
    % --------- 2. Layout Block Mutation ---------
    % بعد از محاسبه numFOPGenes و mutation آنها:
    layoutStart = numFOPGenes + 1;
    offset = 0;
    nMatTotal = 0;

    for t = 1:T
        stageStart = sum([model.Stage(1:t-1).Duration]) + 1;
        stageEnd = stageStart + model.Stage(t).Duration;

        nMat = 0;
        for m = 1:numMaterials
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

            nMatTotal = nMatTotal + 1;
            nMat = nMat + 1;

            % مکان layout این متریال در chromosome
            idx = layoutStart + offset + (nMat - 1) * 5;

            % خواندن FOP جدید
            fop = mutated(nMatTotal);

            % محاسبه ابعاد جدید
            [lx, ly] = GetFootprintFromFOP(model.materials.(mname), fop, relevantDemand, t, model);

            % توجه به زاویه
            theta = mutated(idx + 2);
            if mod(theta, 180) == 90
                [lx, ly] = deal(ly, lx);
            end

            % به‌روزرسانی در layout block
            mutated(idx + 3) = lx;
            mutated(idx + 4) = ly;
        end

        % شمارش Facilityهای فعال برای offset
        nFac = 0;
        facNames = facNames(startsWith(facNames, 'F'));  % فقط آنهایی که با F شروع می‌کنند
        for f = 1:length(facNames)
            fname = facNames{f};
            fobj = model.Facilities.(fname);
            if isfield(fobj, 'activestages') && ismember(t, fobj.activestages)
                nFac = nFac + 1;
            end
        end

        % به‌روزرسانی offset برای stage بعدی
        offset = offset + 5 * (nMat + nFac);
            % --------- 2. Layout Mutation ---------
        M=randi([1 3]);

        switch M
            case 1
                % Swap
                mutated=DoSwap(mutated, model, numFOPGenes, t);

            case 2
                % Reversion
                mutated=DoReversion(mutated, model, numFOPGenes, t);

            case 3
                % Insertion
                mutated=DoInsertion(mutated, model, numFOPGenes, t);
        end
    end
    % --- Repair Step ---
    [mutated, feasible] = RepairChild(mutated, model);
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
function mutated = DoSwap(mutated, model, numFOPGenes, t)

    matNames = fieldnames(model.materials);
    facNames = fieldnames(model.Facilities);
    pos = GetLayoutStartPos(model, t, numFOPGenes);
    
        % --- مرحله 1: لیست مواد فعال در این استیج ---
        stageStart = sum([model.Stage(1:t-1).Duration]) + 1;
        stageEnd = stageStart + model.Stage(t).Duration;

        activeMats = {};
        for m = 1:length(matNames)
            mname = matNames{m};
            demand = model.materials.(mname).DemandFlow;
            if any(demand(stageStart:min(stageEnd, end)))
                activeMats{end+1} = mname;
            end
        end

        % --- مرحله 2: لیست فسیلیتی‌های فعال در این استیج ---
        activeFacs = {};
        for f = 1:length(facNames)
            fname = facNames{f};
            fac = model.Facilities.(fname);
            if ismember(t, fac.activestages) && startsWith(fname, 'F')
                activeFacs{end+1} = fname;
            end
        end

        % --- مرحله 3: لیست ترکیبی آیتم‌ها ---
        allItems = [activeMats, activeFacs];
        nItems = length(allItems);

        % --- مرحله 4: انتخاب تصادفی دو آیتم مختلف ---
        maxTries = 100;
        for attempt = 1:maxTries
            idxs = randsample(nItems, 2);
            idx1 = idxs(1);
            idx2 = idxs(2);

            if ~IsFixedFacility(idx1, t, model) && ~IsFixedFacility(idx2, t, model)
                break;  % موفق شدیم
            end
        end
        % --- مرحله 5: swap در layout block ---
        pos1 = pos + (idx1 - 1) * 5;
        pos2 = pos + (idx2 - 1) * 5;

        tempX = mutated(pos1);
        tempY = mutated(pos1 + 1);
        mutated(pos1)     = mutated(pos2);
        mutated(pos1 + 1) = mutated(pos2 + 1);
        mutated(pos2)     = tempX;
        mutated(pos2 + 1) = tempY;
end
function mutated = DoReversion(mutated, model, numFOPGenes, t)

    matNames = fieldnames(model.materials);
    facNames = fieldnames(model.Facilities);
    pos = GetLayoutStartPos(model, t, numFOPGenes);

        % --- مرحله 1: مواد فعال در این استیج ---
        stageStart = sum([model.Stage(1:t-1).Duration]) + 1;
        stageEnd = stageStart + model.Stage(t).Duration;

        activeMats = {};
        for m = 1:length(matNames)
            mname = matNames{m};
            demand = model.materials.(mname).DemandFlow;
            if any(demand(stageStart:min(stageEnd, end)))
                activeMats{end+1} = mname;
            end
        end

        % --- مرحله 2: فسیلیتی‌های فعال در این استیج ---
        activeFacs = {};
        for f = 1:length(facNames)
            fname = facNames{f};
            fac = model.Facilities.(fname);
            if ismember(t, fac.activestages) && startsWith(fname, 'F')
                activeFacs{end+1} = fname;
            end
        end

        % --- مرحله 3: لیست کامل آیتم‌ها ---
        allItems = [activeMats, activeFacs];
        nItems = length(allItems);

        if nItems < 2
            pos = pos + 5 * nItems;
            % continue;
        end

        % --- مرحله 4: انتخاب بازه تصادفی ---
        maxTries = 100;
        for attempt = 1:maxTries
            idxs = sort(randsample(nItems, 2));  % [i1, i2]
            i1 = idxs(1);
            i2 = idxs(2);

            % بررسی کل بازه [i1:i2] برای وجود Fixed Facility
            hasFixed = false;
            for j = i1:i2
                if IsFixedFacility(j, t, model)
                    hasFixed = true;
                    break;
                end
            end

            if ~hasFixed
                break;  % بازه مناسبه
            end
            if attempt == maxTries
                mutated = DoSwap(mutated, model, numFOPGenes, t);  % fallback
                return;
            end
        end
        % --- مرحله 5: معکوس کردن آیتم‌های بازه [i1:i2] ---
        for k = 0:(i2 - i1)
            src1 = pos + (i1 + k - 1) * 5;
            src2 = pos + (i2 - k - 1) * 5;

            if src1 >= src2
                break;
            end

            tmpX = mutated(src1);
            tmpY = mutated(src1 + 1);
            mutated(src1)     = mutated(src2);
            mutated(src1 + 1) = mutated(src2 + 1);
            mutated(src2)     = tmpX;
            mutated(src2 + 1) = tmpY;
        end
end
function mutated = DoInsertion(mutated, model, numFOPGenes, t)


    matNames = fieldnames(model.materials);
    facNames = fieldnames(model.Facilities);
    pos = GetLayoutStartPos(model, t, numFOPGenes);

        % --- لیست مواد فعال ---
        stageStart = sum([model.Stage(1:t-1).Duration]) + 1;
        stageEnd = stageStart + model.Stage(t).Duration;

        activeMats = {};
        for m = 1:length(matNames)
            mname = matNames{m};
            demand = model.materials.(mname).DemandFlow;
            if any(demand(stageStart:min(stageEnd, end)))
                activeMats{end+1} = mname;
            end
        end

        % --- لیست فسیلیتی‌های فعال ---
        activeFacs = {};
        for f = 1:length(facNames)
            fname = facNames{f};
            fac = model.Facilities.(fname);
            if ismember(t, fac.activestages) && startsWith(fname, 'F')
                activeFacs{end+1} = fname;
            end
        end

        allItems = [activeMats, activeFacs];
        nItems = length(allItems);

        if nItems < 2
            pos = pos + 5 * nItems;
            % continue;
        end

        % --- انتخاب i1 و i2 ---
        maxTries = 100;
        for attempt = 1:maxTries
            idxs = sort(randsample(nItems, 2));  % [i1, i2]
            i1 = idxs(1);
            i2 = idxs(2);

            % بررسی کل بازه [i1:i2] برای وجود Fixed Facility
            hasFixed = false;
            for j = i1:i2
                if IsFixedFacility(j, t, model)
                    hasFixed = true;
                    break;
                end
            end

            if ~hasFixed
                break;  % بازه مناسبه
            end
            if attempt == maxTries
                mutated = DoSwap(mutated, model, numFOPGenes, t);  % fallback
                return;
            end
        end
        % استخراج کل x و y
        layoutXY = zeros(nItems, 2);
        for i = 1:nItems
            layoutXY(i, :) = mutated(pos + (i-1)*5 : pos + (i-1)*5 + 1);
        end
    
        item = layoutXY(i1, :);
        layoutXY(i1, :) = [];
    
        if i1 < i2
            layoutXY = [layoutXY(1:i2-1, :); item; layoutXY(i2:end, :)];
        else
            layoutXY = [layoutXY(1:i2, :); item; layoutXY(i2+1:end, :)];
        end
    
        % بازنویسی فقط x و y
        for i = 1:nItems
            mutated(pos + (i-1)*5)     = layoutXY(i, 1);  % x
            mutated(pos + (i-1)*5 + 1) = layoutXY(i, 2);  % y
        end
end

function pos = GetLayoutStartPos(model, t, numFOPGenes)
    matNames = fieldnames(model.materials);
    facNames = fieldnames(model.Facilities);
    pos = numFOPGenes + 1;  % از پایان FOPها شروع می‌کنیم
    
    for k = 1:t-1
        nMat = 0;
        stageStart = sum([model.Stage(1:k-1).Duration]) + 1;
        stageEnd = stageStart + model.Stage(k).Duration;

        for m = 1:length(matNames)
            demand = model.materials.(matNames{m}).DemandFlow;
            if any(demand(stageStart:min(stageEnd, end)))
                nMat = nMat + 1;
            end
        end

        nFac = 0;
        for f = 1:length(facNames)
            fname = facNames{f};
            if startsWith(fname, 'F') && ismember(k, model.Facilities.(fname).activestages)
                nFac = nFac + 1;
            end
        end

        pos = pos + 5 * (nMat + nFac);
    end
end

function isFixed = IsFixedFacility(idx, t, model)

    matNames = fieldnames(model.materials);
    facNames = fieldnames(model.Facilities);

    stageStart = sum([model.Stage(1:t-1).Duration]) + 1;
    stageEnd = stageStart + model.Stage(t).Duration;

    % --- لیست مواد فعال ---
    activeMats = {};
    for m = 1:length(matNames)
        mname = matNames{m};
        demand = model.materials.(mname).DemandFlow;
        if any(demand(stageStart:min(stageEnd, end)))
            activeMats{end+1} = mname;
        end
    end

    % --- لیست فسیلیتی‌های فعال ---
    activeFacs = {};
    for f = 1:length(facNames)
        fname = facNames{f};
        if startsWith(fname, 'F')
            fac = model.Facilities.(fname);
            if ismember(t, fac.activestages)
                activeFacs{end+1} = fname;
            end
        end
    end

    allItems = [activeMats, activeFacs];

    % اگر ایندکس خارج از بازه بود، false
    if idx < 1 || idx > numel(allItems)
        isFixed = false;
        return;
    end

    name = allItems{idx};

    if startsWith(name, 'M') || startsWith(name, 'B')  % فقط facilities رو چک می‌کنیم
        isFixed = false;
        return;
    end

    % چک کردن Type
    if isfield(model.Facilities.(name), 'Type') && strcmpi(model.Facilities.(name).Type, 'S')
        isFixed = true;
    else
        isFixed = false;
    end
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

        % مرحله 4: قرار دادن B ها روی سایت
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


        % مرحله 5: فسیلیتی‌های ثابت
        if isfield(sol.Stage(t), 'Facilities')
            for i = 1:length(facNames)
                fname = facNames{i};
                if isfield(model.Facilities.(fname), 'Type') && strcmpi(model.Facilities.(fname).Type, 'S')
                    f = sol.Stage(t).Facilities.(fname);
                    fx = max(1, floor(f.x));
                    fy = max(1, floor(f.y));
                    fw = floor(f.Lx);
                    fh = floor(f.Ly);
                    gridMap(fy:fy+fh-1, fx:fx+fw-1) = 1;
                    blockedZones(end+1) = struct('x', fx, 'y', fy, 'width', fw, 'height', fh);
                end
            end
        end

        % مرحله 6: فسیلیتی‌های متحرک
        if isfield(sol.Stage(t), 'Facilities')
            for i = 1:length(facNames)
                fname = facNames{i};
                if isfield(model.Facilities.(fname), 'Type') && strcmpi(model.Facilities.(fname).Type, 'M')
                    f = sol.Stage(t).Facilities.(fname);
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
        end

        % مرحله 7: مواد
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
                                    m.x = x;
                                    m.y = y;
                                    m.Lx = Lx;
                                    m.Ly = Ly;
                                    m.FOP = fop;
                                    sol.Stage(t).Materials.(mname) = m;
                                    placed_this_material = true;
                                    break;
                                else
                                    [x, y] = randomWalk(x, y, gridW, gridH, Lx, Ly);
                                end
                            end
                            if placed_this_material
                                break;
                            end
                        end
                    end
                    if ~placed_this_material
                        allPlaced = false;
                        break;
                    end
                end
                if allPlaced
                    Placed = true;
                else
                    feasible = false;
                    break;
                end
            end
        end
    end


    % === مرحله نهایی: بازسازی کروموزوم
    sol.pos = [];

    % --- FOP
    for t = 1:T
        if isfield(sol.Stage(t), 'Materials')
            matNames = fieldnames(sol.Stage(t).Materials);
            for m = 1:length(matNames)
                mname = matNames{m};
                mat = sol.Stage(t).Materials.(mname);
                if isfield(mat, 'FOP') && ~isempty(mat.FOP)
                    sol.pos(end+1) = mat.FOP;
                end
            end
        end
    end

    % --- Layout: مواد و فسیلیتی‌ها
    for t = 1:T
        if isfield(sol.Stage(t), 'Materials')
            matNames = fieldnames(sol.Stage(t).Materials);
            for m = 1:length(matNames)
                mname = matNames{m};
                mat = sol.Stage(t).Materials.(mname);
                sol.pos(end+1:end+5) = [mat.x, mat.y, mat.theta, mat.Lx, mat.Ly];
            end
        end

        if isfield(sol.Stage(t), 'Facilities')
            facNames = fieldnames(sol.Stage(t).Facilities);
            fFiltered = facNames(startsWith(facNames, 'F'));
            fNumbers = cellfun(@(s) sscanf(s, 'F%d'), fFiltered);
            [~, sortIdx] = sort(fNumbers);
            sortedF = fFiltered(sortIdx);
            for j = 1:length(sortedF)
                fname = sortedF{j};
                f = sol.Stage(t).Facilities.(fname);
                sol.pos(end+1:end+5) = [f.x, f.y, f.theta, f.Lx, f.Ly];
            end
        end
    end

    child_fixed = sol.pos;
end
