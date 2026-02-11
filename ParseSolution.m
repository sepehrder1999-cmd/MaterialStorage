function sol = ParseSolution(pos, model)

    sol.pos = pos;  % ذخیره موقعیت ژن‌ها در ساختار راه‌حل
    % sol.Stage = model.Stage;  % کپی ساختار اولیه stageها

    matNames = fieldnames(model.materials);
    facNames = fieldnames(model.Facilities);
    T = numel(model.Stage);
    numMaterials = length(matNames);

    k = 1;

    for i=1:T    
        sol.Stage(i).occupiedRects = struct( ...
                'name', {}, ...
                'x', {}, ...
                'y', {}, ...
                'width', {}, ...
                'height', {});

        facilityNames = fieldnames(model.Facilities);
        activeFacs = {};
        for j = 1:length(facilityNames)
            fname = facilityNames{j};
            f = model.Facilities.(fname);
            if ismember(i, f.activestages)
                activeFacs{end+1} = fname;
            end
        end


        % ساختمان‌های ثابت (B1، B2، G)
        fname = startsWith(facilityNames, {'B', 'G'});
        fname = facilityNames(fname);  % فقط B و Gها
        for j = 1:numel(fname)
            name = fname{j};  % مقدار رشته‌ای از cell array
                
            f = model.Facilities.(name);
                
            if ismember(i, f.activestages)
                x = f.x; y = f.y; w = f.Lx; h = f.Ly;
                sol.Stage(i).Buildings.(name) = struct('x',x,'y',y,'Lx',w,'Ly',h);
                 sol.Stage(i).occupiedRects(end+1) = struct( ...
                    'name', name, ...
                    'x', x, ...
                    'y', y, ...
                    'width', w, ...
                    'height', h);

            end
        end
    end
    
    % --------- 1. Parse FOP Block ---------
    for t = 1:T
        for m = 1:numMaterials
            mname = matNames{m};
            demand = model.materials.(mname).DemandFlow;
    
            stageStart = sum([model.Stage(1:t-1).Duration]) + 1;
            stageEnd = stageStart + model.Stage(t).Duration;
    
            if any(demand(stageStart:min(stageEnd, end)))
                % --- گرفتن FOP از ژنوم ---
                fop = pos(k);
                sol.Stage(t).Materials.(mname).FOP = fop;
                k = k + 1;
    
                % --- محاسبه orderQty ---
                relevantDemand = demand(stageStart:min(stageEnd, end));
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
                            orderQty = [relevantDemand(2)-demand(stageStart); diff(relevantDemand)];
                        else
                            orderQty = [];
                        end
                    else
                        for s = 1:nSteps
                            w1 = (s-1)*fop + 1;
                            w2 = s*fop + 1;
        
                            if w1 > length(relevantDemand)
                                continue; % این بازه دیگه در دسترس نیست
                            end
        
                            if w2 > length(relevantDemand)
                                w2 = length(relevantDemand);
                            end
        
                            % حالت خاص: اولین سفارش در اولین استیج
                            if t == 1 && s == 1
                                orderQty(s) = relevantDemand(w2);
                                continue;
                            end
        
                            orderQty(s) = relevantDemand(w2) - relevantDemand(w1);
                        end
                    end
                end
    
                % --- ذخیره در sol ---
                sol.Stage(t).Materials.(mname).orderQty = orderQty;
            end
        end
    end

    % --------- 2. Parse Layout Block ---------
    for t = 1:T
        % --- Materials ---
        for m = 1:numMaterials
            mname = matNames{m};
            demand = model.materials.(mname).DemandFlow;

            stageStart = sum([model.Stage(1:t-1).Duration]) + 1;
            stageEnd = stageStart + model.Stage(t).Duration;

            if any(demand(stageStart:min(stageEnd, end)))
                sol.Stage(t).Materials.(mname).x     = pos(k);
                sol.Stage(t).Materials.(mname).y     = pos(k+1);
                sol.Stage(t).Materials.(mname).theta = pos(k+2);
                sol.Stage(t).Materials.(mname).Lx    = pos(k+3);
                sol.Stage(t).Materials.(mname).Ly    = pos(k+4);
                sol.Stage(t).occupiedRects(end+1) = struct( ...
                'name', mname, ...
                'x', sol.Stage(t).Materials.(mname).x, ...
                'y', sol.Stage(t).Materials.(mname).y, ...
                'width', sol.Stage(t).Materials.(mname).Lx, ...
                'height', sol.Stage(t).Materials.(mname).Ly);

                k = k + 5;
            end
        end

        % --- Facilities ---
        % فقط آن‌هایی که با F شروع می‌شن
        filteredF = facNames(startsWith(facNames, 'F'));

        % مرتب‌سازی بر اساس شماره
        fNumbers = cellfun(@(s) sscanf(s, 'F%d'), filteredF);
        [~, sortIdx] = sort(fNumbers);
        sortedF = filteredF(sortIdx);

        for j = 1:length(sortedF)
            fname = sortedF{j};
            f = model.Facilities.(fname);
            if ismember(t, f.activestages)
                sol.Stage(t).Facilities.(fname).x     = pos(k);
                sol.Stage(t).Facilities.(fname).y     = pos(k+1);
                sol.Stage(t).Facilities.(fname).theta = pos(k+2);
                sol.Stage(t).Facilities.(fname).Lx    = pos(k+3);
                sol.Stage(t).Facilities.(fname).Ly    = pos(k+4);
                sol.Stage(t).occupiedRects(end+1) = struct( ...
                'name', fname, ...
                'x', sol.Stage(t).Facilities.(fname).x, ...
                'y', sol.Stage(t).Facilities.(fname).y, ...
                'width', sol.Stage(t).Facilities.(fname).Lx, ...
                'height', sol.Stage(t).Facilities.(fname).Ly);

                k = k + 5;
            end
        end
    end
end