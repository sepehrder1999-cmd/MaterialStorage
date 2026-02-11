function raw_penalties_all_stages = calculateRaw5SPenalties(sol, model)

    sol = ParseSolution(sol,model);
    T = numel(model.Stage); % تعداد کل فازها
    raw_penalties_all_stages = repmat(struct('P_Sort', 0, 'P_SetInOrder', 0, 'P_Shine', 0, 'P_Standardize', 0), T, 1);

    % =================================================================
    %             حلقه اصلی برای محاسبه جریمه‌ها در هر فاز
    % =================================================================
    for t = 1:T
        
        sol_stage = sol.Stage(t); % استخراج اطلاعات فاز فعلی

        % مقداردهی اولیه برای فاز فعلی
        P_Sort = 0;
        P_SetInOrder = 0;
        P_Shine = 0;
        P_Standardize = 0;

        % ======================== محاسبه P_Sort ========================
        P_Sort = 0;
        Ws1 = model.s5_params.ws1;
        Ws2 = model.s5_params.ws2;
        
        if isfield(sol_stage, 'Materials')
            materials = sol_stage.Materials;
            mat_names = fieldnames(materials);
            
            for m = 1:length(mat_names)
                mat_id = mat_names{m};
                mat_storage = materials.(mat_id);
                mat_model = model.materials.(mat_id);
                
                orderQty_vector = mat_storage.orderQty;
                
                if isempty(orderQty_vector)
                    continue;
                end
        
                % --- محاسبه بخش اول: جریمه نگهداری موجودی (صحیح است) ---
                penalty1 = 0;
                deltaT = mat_storage.FOP;
                for i = 1:length(orderQty_vector)
                    qty_at_period_i = orderQty_vector(i);
                    penalty1 = penalty1 + (qty_at_period_i * deltaT);
                end
        
                % --- محاسبه بخش دوم: جریمه فضای اضافه (اصلاح شده) ---
                area_alloc = mat_storage.Lx * mat_storage.Ly;
                
                % پیدا کردن بزرگترین مقدار سفارش
                max_order_qty = max(orderQty_vector);
                
                % محاسبه مساحت مورد نیاز فقط برای بزرگترین سفارش
                area_req_for_max_qty = 0;
                footprintOptions = mat_model.StorageFootprint;
                for f = 1:length(footprintOptions)
                    qrange = footprintOptions(f).Qrange;
                    if max_order_qty >= qrange(1) && max_order_qty <= qrange(2)
                       area_req_for_max_qty = footprintOptions(f).Lx * footprintOptions(f).Ly;
                       break;
                   end
                end
                
                % محاسبه جریمه نهایی فقط یک بار
                if area_req_for_max_qty > 0
                    penalty2 = max(0, area_alloc - area_req_for_max_qty);
                else
                    penalty2 = 0;
                end
                
                % --- تجمیع وزنی دو بخش برای این متریال ---
                P_Sort_material = Ws1 * penalty1 + Ws2 * penalty2;
                                         
                P_Sort = P_Sort + P_Sort_material;
            end
        end
        % ===================== محاسبه P_SetInOrder =====================
        % بخش اول: جریمه فاصله تسهیلات (Adjacency)
        P_Adjacency = 0;
        Wo1 = model.s5_params.wo1;
        Wo2 = model.s5_params.wo2;
        if isfield(sol_stage, 'Facilities')
            facilities = sol_stage.Facilities;
            fac_names = fieldnames(facilities);
            num_facs = length(fac_names);
            
            for i = 1:num_facs
                for j = i + 1:num_facs
                    f_num1 = sscanf(fac_names{i}, 'F%d');
                    f_num2 = sscanf(fac_names{j}, 'F%d');
                    closeness_rating = model.s5_params.ClosenessMatrix(f_num1, f_num2);
                    
                    if closeness_rating > 0
                        fac1 = facilities.(fac_names{i});
                        fac2 = facilities.(fac_names{j});
                        dist = sqrt((fac1.x - fac2.x)^2 + (fac1.y - fac2.y)^2);
                        P_Adjacency = P_Adjacency + (dist * closeness_rating);
                    end
                end
            end
        end
        
        % بخش دوم: جریمه دسترسی به مصالح (Access)
        P_Access = 0;
        if isfield(sol_stage, 'Materials')
            materials = sol_stage.Materials;
            mat_names = fieldnames(materials);
            
            for k = 1:length(mat_names)
                mat_id = mat_names{k}; % e.g., 'M1'
                mat_storage = materials.(mat_id);
                
                % --- یافتن تمام مصرف‌کنندگان و وزن تقاضای آنها ---
                consumer_locations = []; % ماتریس برای ذخیره مختصات X, Y
                consumer_weights = [];   % بردار برای ذخیره وزن‌ها (RequiredQty)
        
                for h_idx = 1:length(model.handlingData)
                    handling_entry = model.handlingData(h_idx);
                    
                    % بررسی تطابق متریال و فاز
                    if strcmp(handling_entry.MaterialID, mat_id) && handling_entry.Stage == t
                        consumer_name = handling_entry.Facility; % e.g., 'B1' or 'F2'
                        
                        % یافتن مختصات مصرف‌کننده در جانمایی فعلی
                        consumer_fac = struct();
                        if isfield(sol_stage.Facilities, consumer_name)
                            consumer_fac = sol_stage.Facilities.(consumer_name);
                        elseif isfield(model.Facilities, consumer_name) % برای ساختمان‌های ثابت B1, B2
                            consumer_fac = model.Facilities.(consumer_name);
                        end
                        
                        if ~isempty(fieldnames(consumer_fac))
                            % محاسبه مرکز تسهیلات مصرف‌کننده
                            center_x = consumer_fac.x + consumer_fac.Lx / 2;
                            center_y = consumer_fac.y + consumer_fac.Ly / 2;
                            
                            % ذخیره مکان و وزن (مقدار تقاضا)
                            consumer_locations = [consumer_locations; center_x, center_y];
                            consumer_weights = [consumer_weights; handling_entry.RequiredQty];
                        end
                    end
                end
                
                % --- محاسبه نقطه مصرف بهینه (مرکز ثقل وزنی) ---
                if ~isempty(consumer_locations)
                    total_weight = sum(consumer_weights);
                    
                    % برای جلوگیری از تقسیم بر صفر اگر مجموع وزن‌ها صفر باشد
                    if total_weight > 1e-6 
                        % محاسبه مختصات X و Y مرکز ثقل وزنی
                        weighted_x = sum(consumer_locations(:,1) .* consumer_weights) / total_weight;
                        weighted_y = sum(consumer_locations(:,2) .* consumer_weights) / total_weight;
                        
                        weighted_consumption_point = [weighted_x, weighted_y];
                    else
                        % اگر وزن‌ها صفر باشند، از مرکز هندسی ساده استفاده کن
                        weighted_consumption_point = mean(consumer_locations, 1);
                    end
        
                    % محاسبه فاصله انبار تا نقطه مصرف بهینه
                    dist = sqrt((mat_storage.x - weighted_consumption_point(1))^2 + ...
                                (mat_storage.y - weighted_consumption_point(2))^2);
                    
                    % در اینجا، جریمه خود فاصله است. وزن‌ها قبلا در محاسبه نقطه مصرف لحاظ شده‌اند.
                    P_Access = P_Access + dist;
                end
            end
        end
        
        P_SetInOrder = Wo1 * P_Adjacency + Wo2 * P_Access; % به‌روزرسانی مقدار نهایی
        % ======================= محاسبه P_Shine =======================
        % --- بخش اول: شاخص تراکم سایت (SCI) ---
        total_occupied_area = 0;
        if isfield(sol_stage, 'Facilities')
            facilities = sol_stage.Facilities;
            fac_names = fieldnames(facilities);
            for i = 1:length(fac_names)
                fac = facilities.(fac_names{i});
                total_occupied_area = total_occupied_area + (fac.Lx * fac.Ly);
            end
        end
        if isfield(sol_stage, 'Materials')
            materials = sol_stage.Materials;
            mat_names = fieldnames(materials);
            for k = 1:length(mat_names)
                mat = materials.(mat_names{k});
                total_occupied_area = total_occupied_area + (mat.Lx * mat.Ly);
            end
        end
        
        Penalty_SCI = total_occupied_area / model.s5_params.UsableArea;
        % --- بخش دوم: جریمه تداخل با مسیرهای تردد ---
        Penalty_PathOverlap = 0;
        corridors = model.s5_params.TrafficCorridors;
        Buildings = fieldnames(sol_stage.Buildings);
        
        allitems = struct();
        is_building = startsWith(Buildings,'B');
        building_names = Buildings(is_building);
        Buildings = sol_stage.Buildings;
        for i = 1:length(building_names)
            bname = building_names{i};
            allitems.(bname) = Buildings.(bname);
        end
        for i = 1:length(fac_names)
            fname = fac_names{i};
            allitems.(fname) = facilities.(fname);
        end
        allitemsnames = fieldnames(allitems);
        for i = 1:length(allitemsnames)
                item = allitems.(allitemsnames{i});
                item_rect = [item.x, item.y, item.Lx, item.Ly];
            for j = 1:size(corridors, 1)
                corridor_rect = corridors(j, :);
                overlap_area = calculateOverlapArea(item_rect, corridor_rect);
                Penalty_PathOverlap = Penalty_PathOverlap + overlap_area;
            end
        end

        % --- بخش سوم: جریمه نقض مناطق ایمنی ---
        Penalty_SafetyViolation = 0;
        safety_zones = model.s5_params.SafetyZones;
        for i = 1:length(safety_zones)
            zone = safety_zones(i);
            zone_rect = [];
            
            % ساخت مستطیل منطقه ایمنی بر اساس تعریف آن
            if isfield(sol_stage.Facilities, zone.center_facility)
                center_fac = sol_stage.Facilities.(zone.center_facility);
                
                if strcmp(zone.type, 'rectangle')
                    zone_rect = [center_fac.x - zone.buffer_x, ...
                                 center_fac.y - zone.buffer_y, ...
                                 center_fac.Lx + 2*zone.buffer_x, ...
                                 center_fac.Ly + 2*zone.buffer_y];
                elseif strcmp(zone.type, 'circle')
                    % تقریب دایره با یک مربع محاطی برای سادگی محاسبه همپوشانی
                    zone_rect = [center_fac.x - zone.radius, ...
                                 center_fac.y - zone.radius, ...
                                 2*zone.radius, ...
                                 2*zone.radius];
                end
            end
            
            if ~isempty(zone_rect)
                % بررسی تداخل این منطقه با تمام آیتم‌های دیگر
                for j = 1:length(allitemsnames)
                    item = allitems.(allitemsnames{j});
                    
                    % یک تسهیلات نباید با منطقه ایمنی خودش تداخل داشته باشد
                    if isfield(item,'Name') && strcmp(item.Name, zone.center_facility)
                        continue;
                    end
        
                    item_rect = [item.x, item.y, item.Lx, item.Ly];
                    overlap_area = calculateOverlapArea(item_rect, zone_rect);
                    Penalty_SafetyViolation = Penalty_SafetyViolation + overlap_area;
                end
            end
        end
        % --- تجمیع وزنی سه بخش برای P_Shine نهایی ---
        P_Shine = model.s5_params.wc1 * Penalty_SCI + ...
          model.s5_params.wc2 * Penalty_PathOverlap + ...
          model.s5_params.wc3 * Penalty_SafetyViolation;
        
        % ==================== محاسبه P_Standardize ====================
        if t > 1
            prev_sol_stage = sol.Stage(t-1);
            if isfield(sol_stage, 'Facilities') && isfield(prev_sol_stage, 'Facilities')
                current_fac_names = fieldnames(sol_stage.Facilities);
                prev_fac_names = fieldnames(prev_sol_stage.Facilities);
                
                % پیدا کردن تسهیلات مشترک بین فاز فعلی و قبلی
                common_fac_names = intersect(current_fac_names, prev_fac_names);
                
                total_distance_moved = 0;
                for i = 1:length(common_fac_names)
                    fname = common_fac_names{i};
                    fac_current = sol_stage.Facilities.(fname);
                    fac_prev = prev_sol_stage.Facilities.(fname);
                    
                    dist_moved = sqrt((fac_current.x - fac_prev.x)^2 + (fac_current.y - fac_prev.y)^2);
                    total_distance_moved = total_distance_moved + dist_moved;
                end
                P_Standardize = total_distance_moved;
            end
        end
        
        % ================== ذخیره نتایج برای فاز فعلی ==================
        raw_penalties_all_stages(t).P_Sort = P_Sort;
        raw_penalties_all_stages(t).P_SetInOrder = P_SetInOrder;
        raw_penalties_all_stages(t).P_Shine = P_Shine;
        raw_penalties_all_stages(t).P_Standardize = P_Standardize;
    end
end
function area = calculateOverlapArea(rect1, rect2)
    % rect format: [x, y, Lx, Ly]
    
    x_overlap = max(0, min(rect1(1) + rect1(3), rect2(1) + rect2(3)) - max(rect1(1), rect2(1)));
    y_overlap = max(0, min(rect1(2) + rect1(4), rect2(2) + rect2(4)) - max(rect1(2), rect2(2)));
    
    area = x_overlap * y_overlap;
end