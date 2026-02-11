function [LC, sol2] = LayoutCost(sol, model, s5_norm_params, solraw)
    
    global NFE;
    if isempty(NFE)
        NFE = 0;
    end
    NFE = NFE + 1;

    T = numel(model.Stage);

    OC = zeros(1, T);  % مقداردهی اولیه برای همه مراحل
    MHC = zeros(1, T);  % مقداردهی اولیه کل MHC‌ها
    RTC = zeros(1, T);  % مقداردهی اولیه برای هر مرحله
    SRC = zeros(1, T);  % مقداردهی اولیه بردار SRC برای تمام مراحل
    penalty = zeros(1, T);  % مقداردهی اولیه برای هر مرحله
    Penalty_5S_per_stage = zeros(1, T);
    S5_Separated = struct('Sort', [], 'SetInOrder', [], 'Shine', [], 'Standardize', [], 'S5_score', []);


    prevFacs = struct();

    sol = ParseSolution(sol, model);


    for i = 1:T
        facPos = struct();
        matPos = struct();
        violatedConstraints = {};  % لاگ نقض‌ها

        
        facNames = [ fieldnames(sol.Stage(i).Facilities); fieldnames(sol.Stage(i).Buildings) ];
        for f = 1:numel(facNames)
            fname = facNames{f};
                % تشخیص نوع منبع داده
            if startsWith(fname, 'F')
                fdata = sol.Stage(i).Facilities.(fname);
            else 
                fdata = sol.Stage(i).Buildings.(fname);
            end
            cx = fdata.x + fdata.Lx/2;
            cy = fdata.y + fdata.Ly/2;

            if isfield(fdata, 'theta')
                    theta = fdata.theta;
                else
                    theta = 0;   % پیش‌فرض
            end
            facPos.(fname) = [cx, cy, theta];
        end

        if isfield(sol.Stage(i), 'Materials')
            matNames = fieldnames(sol.Stage(i).Materials);
            for m = 1:numel(matNames)
                mname = matNames{m};
                mdata = sol.Stage(i).Materials.(mname);
                cx = mdata.x + mdata.Lx/2;
                cy = mdata.y + mdata.Ly/2;
                matPos.(mname) = [cx, cy];
            end
        end

    %% Calculating Ordering Cost (OC)
    
        for t = 1:T
            stageOC = 0;
    
            if isfield(sol.Stage(t), 'Materials')
                matNames = fieldnames(sol.Stage(t).Materials);
    
                for m = 1:numel(matNames)
                    mname = matNames{m};
                    mdata = sol.Stage(t).Materials.(mname);
    
                    if ~isfield(mdata, 'orderQty')
                        continue;   % اگر هنوز orderQty ست نشده
                    end
    
                    orders = mdata.orderQty;  % می‌تونه یه عدد یا آرایه باشه
    
                    for q = orders                        
                        PCR = 0;
                        purchaseTable = model.materials.(mname).PurchaseCost;
                        for k = 1:numel(purchaseTable)
                            row = purchaseTable(k);
                            if any(row.Stage == t) && q >= row.Quantity(1) && q <= row.Quantity(2)
                                PCR = row.PCR;
                                break;
                            end
                        end
    
                        % --- پیدا کردن DLC بر اساس بازه Q ---
                        DLC = 0;
                        deliveryTable = model.materials.(mname).DeliveryCost;
                        for k = 1:numel(deliveryTable)
                            row = deliveryTable(k);
                            if q >= row.Quantity(1) && q <= row.Quantity(2)
                                DLC = row.DLC;
                                break;
                            end
                        end
    
                        % --- هزینه سفارش ---
                        stageOC = stageOC + (q * PCR + DLC);
                    end
                end
            end
    
            OC(t) = stageOC;
        end

        % --- ساخت یک استراکت برای ذخیره خروجی ---
        % ساختار خروجی: آرایه‌ای از استراکت‌ها
        OrderQty = struct('Stage', {}, 'MatID', {}, 'FOP', {}, 'OrderQty', {});
        entryIndex = 1;
        
        % --- حلقه اصلی برای گردش روی تمام استیج‌ها ---
        for t = 1:T
            % بررسی اینکه آیا متریالی در استیج فعلی وجود دارد یا خیر
            if isfield(sol.Stage(t), 'Materials')
                
                % استخراج نام تمام متریال‌های موجود در استیج t
                matNames = fieldnames(sol.Stage(t).Materials);
        
                % --- حلقه داخلی برای گردش روی هر متریال ---
                for m = 1:numel(matNames)
                    mname = matNames{m};
                    mdata = sol.Stage(t).Materials.(mname);
        
                    % بررسی وجود فیلدهای مورد نظر برای جلوگیری از خطا
                    if isfield(mdata, 'FOP') && isfield(mdata, 'orderQty')
                        
                        % استخراج مقادیر FOP و orderQty
                        fopValue = mdata.FOP;
                        orderQuantities = mdata.orderQty;
        
                        % --- ذخیره اطلاعات در یک عنصر جدید از آرایه استراکت ---
                        OrderQty(entryIndex).Stage = t;
                        OrderQty(entryIndex).MatID = mname;
                        OrderQty(entryIndex).FOP = fopValue;
                        OrderQty(entryIndex).OrderQty = orderQuantities;
                        
                        entryIndex = entryIndex + 1;
                    end
                end
            end
        end        

    %% Calculating MHC
        for c = 1:length(model.handlingData)
            Q = model.handlingData(c).RequiredQty;
            q = model.handlingData(c).HandlingQty;
            HCR = model.handlingData(c).HourlyCost;
            v = model.handlingData(c).Speed;
        
            C = (2 * (Q / q) * HCR) / v;
            model.handlingData(c).TravelCostRate = C;
        end

        stageRows = find([model.handlingData.Stage] == i);
        stageCost = 0;
        for idx = stageRows
            row = model.handlingData(idx);
            matID = row.MaterialID;
            facName = row.Facility;
        
            if ~isfield(matPos, matID) || ~isfield(facPos, facName)
                warning('Missing position for %s or %s', matID, facName);
                continue;
            end
        
            dist = norm(matPos.(matID) - facPos.(facName)(1:2));
            rowCost = row.TravelCostRate * dist;
            stageCost = stageCost + rowCost;
        end
        MHC(i) = stageCost;

    %% Calculating RTC
        activeNow = facNames; % اینجا میشه customize کرد
        stageRTC = 0;
        for a = 1:numel(activeNow)
            for b = a+1:numel(activeNow)
                fa = activeNow{a};
                fb = activeNow{b};
                if ~isfield(facPos, fa) || ~isfield(facPos, fb)
                    continue;
                end
                dist = norm(facPos.(fa)(1:2) - facPos.(fb)(1:2));
            
                idA = find(strcmp(model.facilityIndex, fa));
                idB = find(strcmp(model.facilityIndex, fb));
                if isempty(idA) || isempty(idB), continue; end
                    C = model.travelCostRate(idA, idB);
                    if isnan(C) || C == 0, continue; end
                    stageRTC = stageRTC + C * dist;
             end
         end
        RTC(i) = stageRTC;

    %% Calculating SRC
        if i == 1
            SRC(i) = 0;
        else
            SRC_t = 0;
            facNamesNow = fieldnames(facPos);
            for f = 1:numel(facNamesNow)
                fname = facNamesNow{f};
                if isfield(prevFacs, fname)
                    dist = norm(facPos.(fname)(1:2) - prevFacs.(fname)(1:2));
                    thetaChanged = facPos.(fname)(3) ~= prevFacs.(fname)(3);
                    if (dist > 0 || thetaChanged) && isfield(model.Facilities.(fname), 'RelCost')
                        SRC_t = SRC_t + model.Facilities.(fname).RelCost;
                    end
                end
            end
        SRC(i) = SRC_t;
        end
    
        prevFacs = facPos;

    %% Penalty for Constraints

        Constraints = model.Constraints;
        thisPenalty = 0;  % برای مرحله i
        
    
        % Safety (Min Distance)
        for s = 1:size(Constraints.safetyMinPairs,1)
            A = Constraints.safetyMinPairs{s,1};
            B = Constraints.safetyMinPairs{s,2};
            dmin = Constraints.safetyMinPairs{s,3};
        
            % گرفتن داده کامل آیتم (نه فقط مرکز)
            if isfield(sol.Stage(i).Facilities, A)
                Adata = sol.Stage(i).Facilities.(A);
            elseif isfield(sol.Stage(i).Buildings, A)
                Adata = sol.Stage(i).Buildings.(A);
            elseif isfield(sol.Stage(i).Materials, A)
                Adata = sol.Stage(i).Materials.(A);
            else
                continue;
            end        

            if isfield(sol.Stage(i).Facilities, B)
                Bdata = sol.Stage(i).Facilities.(B);
            elseif isfield(sol.Stage(i).Buildings, B)
                Bdata = sol.Stage(i).Buildings.(B);
            elseif isfield(sol.Stage(i).Materials, B)
                Bdata = sol.Stage(i).Materials.(B);
            else
                continue;
            end       

            dist = rectDistance(Adata, Bdata);
            if dist < dmin
                if (dmin - dist) <1
                    thisPenalty = thisPenalty + 50000 * (dmin - dist);
                elseif (dmin - dist) <3
                    thisPenalty = thisPenalty + 100000 * (dmin - dist);
                else
                    thisPenalty = thisPenalty + 200000 * (dmin - dist);
                end
                violatedConstraints{end+1} = sprintf('SafetyMin violated: %s-%s (dist=%.2f < %.2f)', ...
                                             A, B, dist, dmin);
            end
        end
        
        % Operational (Max Distance)
        for s = 1:size(Constraints.operationalMaxPairs,1)
            A = Constraints.operationalMaxPairs{s,1};
            B = Constraints.operationalMaxPairs{s,2};
            dmax = Constraints.operationalMaxPairs{s,3};
        
            if isfield(facPos, A)
                posA = facPos.(A);
            elseif isfield(matPos, A)
                posA = matPos.(A);
            else
                continue;
            end
        
            if isfield(facPos, B)
                posB = facPos.(B);
            elseif isfield(matPos, B)
                posB = matPos.(B);
            else
                continue;
            end
        
            dist = norm(posA(1:2) - posB(1:2));
            if dist > dmax
                if (dist - dmax) <1
                    thisPenalty = thisPenalty + 50000 * (dist - dmax);
                elseif (dist - dmax) <3
                    thisPenalty = thisPenalty + 100000 * (dist - dmax);
                else
                    thisPenalty = thisPenalty + 200000 * (dist - dmax);
                end
                violatedConstraints{end+1} = sprintf('OperationalMax violated: %s-%s (dist=%.2f > %.2f)', ...
                                             A, B, dist, dmax);
            end
        end
        
        % Operational (Min Distance)
        for s = 1:size(Constraints.operationalMinPairs,1)
            A = Constraints.operationalMinPairs{s,1};
            B = Constraints.operationalMinPairs{s,2};
            dmin = Constraints.operationalMinPairs{s,3};
        
            if isfield(sol.Stage(i).Facilities, A)
                Adata = sol.Stage(i).Facilities.(A);
            elseif isfield(sol.Stage(i).Buildings, A)
                Adata = sol.Stage(i).Buildings.(A);
            elseif isfield(sol.Stage(i).Materials, A)
                Adata = sol.Stage(i).Materials.(A);
            else
                continue;
            end
        
            if isfield(sol.Stage(i).Facilities, B)
                Bdata = sol.Stage(i).Facilities.(B);
            elseif isfield(sol.Stage(i).Buildings, B)
                Bdata = sol.Stage(i).Buildings.(B);
            elseif isfield(sol.Stage(i).Materials, B)
                Bdata = sol.Stage(i).Materials.(B);
            else
                continue;
            end        
            dist = rectDistance(Adata, Bdata);    % فاصله مرزی
            if dist < dmin
                if (dmin - dist) <1
                    thisPenalty = thisPenalty + 50000 * (dmin - dist);
                elseif (dmin - dist) <3
                    thisPenalty = thisPenalty + 100000 * (dmin - dist);
                else
                    thisPenalty = thisPenalty + 200000 * (dmin - dist);
                end
                violatedConstraints{end+1} = sprintf('OperationalMin violated: %s-%s (dist=%.2f < %.2f)', ...
                                             A, B, dist, dmin);
            end
        end
        
        penalty(i) = thisPenalty;
        sol2.violations{i} = violatedConstraints;

        % ===================================================================
        %         بخش اضافه شده برای محاسبه هزینه نهایی با جریمه 5S
        % ===================================================================
        % دریافت مقادیر خام که در GA.m محاسبه و به ساختار sol اضافه شد
        % فقط برای فاز اول (برای سادگی)
        raw_penalties = solraw.raw_penalties(i); 

            % --- نرمال‌سازی (رابطه 7) ---
        % برای جلوگیری از تقسیم بر صفر
        if (s5_norm_params.Sort_max - s5_norm_params.Sort_min) > 1e-6
            P_Sort_norm = (raw_penalties.P_Sort - s5_norm_params.Sort_min) / ...
                           (s5_norm_params.Sort_max - s5_norm_params.Sort_min);
        else
            P_Sort_norm = 0;
        end

        if (s5_norm_params.SetInOrder_max - s5_norm_params.SetInOrder_min) > 1e-6
            P_SetInOrder_norm = (raw_penalties.P_SetInOrder - s5_norm_params.SetInOrder_min) / ...
                        (s5_norm_params.SetInOrder_max - s5_norm_params.SetInOrder_min);
        else
            P_SetInOrder_norm = 0;
        end

        if (s5_norm_params.Shine_max - s5_norm_params.Shine_min) > 1e-6
            P_Shine_norm = (raw_penalties.P_Shine - s5_norm_params.Shine_min) / ...
                       (s5_norm_params.Shine_max - s5_norm_params.Shine_min);
        else
            P_Shine_norm = 0;
        end

        if (s5_norm_params.Standardize_max - s5_norm_params.Standardize_min) > 1e-6
            P_Standardize_norm = (raw_penalties.P_Standardize - s5_norm_params.Standardize_min) / ...
                       (s5_norm_params.Standardize_max - s5_norm_params.Standardize_min);
        else
            P_Standardize_norm = 0;
        end

            % --- تجمیع وزنی (رابطه 8) ---
            Penalty_5S = model.s5_params.Ws * P_Sort_norm + ...
                 model.s5_params.Wo * P_SetInOrder_norm + ...
                 model.s5_params.Wc * P_Shine_norm + ...
                 model.s5_params.Wt * P_Standardize_norm;

            % --- کد اضافه شده: ذخیره جریمه 5S برای فاز i ---
            Penalty_5S_per_stage(i) = Penalty_5S;

            S5_Separated(i).Sort = P_Sort_norm;
            S5_Separated(i).SetInOrder = P_SetInOrder_norm;
            S5_Separated(i).Shine = P_Shine_norm;
            S5_Separated(i).Standardize = P_Standardize_norm;
            S5_Separated(i).S5_score = Penalty_5S_per_stage(i);




    

    end

    CLC = OC + MHC + RTC + SRC + penalty;

    % --- تابع هدف نهایی با ضریب آلفا ---
    alpha = model.s5_params.alpha;
    Z = CLC + alpha * Penalty_5S_per_stage;
    LC = sum(Z(:));
    
    % ===================================================================
    
        
    sol2.OC  = OC;
    sol2.MHC=MHC;
    sol2.RTC=RTC;
    sol2.SRC=SRC;
    sol2.OrderQty=OrderQty;
    sol2.penalty = penalty;
    sol2.S5_Separated = S5_Separated;
    sol2.Penalty_5S = alpha * Penalty_5S_per_stage;
    sol2.Z=Z;
    sol2.LC=LC;


end
function d = rectDistance(A, B)
    % A و B ساختار با فیلدهای x, y, Lx, Ly هستند
    dx = max(0, max(B.x - (A.x + A.Lx), A.x - (B.x + B.Lx)));
    dy = max(0, max(B.y - (A.y + A.Ly), A.y - (B.y + B.Ly)));
    d = sqrt(dx^2 + dy^2);
end
