function sol = ConstraintRepair(sol, model, it)

    if it < 40
        return;
    end

    constraints = model.Constraints;
    stages = sol.Stage;

    for i = 1:numel(stages)
        rects = stages(i).occupiedRects;

        % --- Safety (Min Distance) ---
        for s = 1:size(constraints.safetyMinPairs,1)
            A = constraints.safetyMinPairs{s,1};
            B = constraints.safetyMinPairs{s,2};
            dmin = constraints.safetyMinPairs{s,3};

            [rectA, rectB] = getRects(rects, A, B);
            if isempty(rectA) || isempty(rectB), continue; end

            d = rectDistance(rectA, rectB);
            if d < dmin
                sol = tryRelocate(sol, model, i, A, rectA);
                sol = tryRelocate(sol, model, i, B, rectB);
            end
        end

        % --- Operational (Max Distance) ---
        for s = 1:size(constraints.operationalMaxPairs,1)
            A = constraints.operationalMaxPairs{s,1};
            B = constraints.operationalMaxPairs{s,2};
            dmax = constraints.operationalMaxPairs{s,3};

            [rectA, rectB] = getRects(rects, A, B);
            if isempty(rectA) || isempty(rectB), continue; end

            % فاصله مرکز به مرکز
            d = norm([rectA.x + rectA.width/2, rectA.y + rectA.height/2] - ...
                  [rectB.x + rectB.width/2, rectB.y + rectB.height/2]);
            if d > dmax
                sol = tryRelocate(sol, model, i, A, rectA);
                sol = tryRelocate(sol, model, i, B, rectB);
            end
        end

        % --- Operational (Min Distance) ---
        for s = 1:size(constraints.operationalMinPairs,1)
            A = constraints.operationalMinPairs{s,1};
            B = constraints.operationalMinPairs{s,2};
            dmin = constraints.operationalMinPairs{s,3};

            [rectA, rectB] = getRects(rects, A, B);
            if isempty(rectA) || isempty(rectB), continue; end

            d = rectDistance(rectA, rectB);
            if d < dmin
                sol = tryRelocate(sol, model, i, A, rectA);
                sol = tryRelocate(sol, model, i, B, rectB);
            end
        end
    end
end

%% تابع کمکی: پیدا کردن مستطیل از اسم
function [rectA, rectB] = getRects(rects, A, B)
    rectA = [];
    rectB = [];
    for r = 1:numel(rects)
        if ~isfield(rects(r), 'name')
            continue; % اگر name وجود نداشت، رد شو
        end
        thisName = rects(r).name;
        if isstring(thisName) || ischar(thisName)
            thisName = char(thisName); % همیشه به رشته تبدیل کن
        end

        if strcmp(thisName, A)
            rectA = rects(r);
        elseif strcmp(thisName, B)
            rectB = rects(r);
        end
    end
end

%% تابع کمکی: اگر قابل جابجایی بود → randomWalk
function sol = tryRelocate(sol, model, stageIdx, name, rect)

    if isempty(rect), return; end

    % Facilities
    if isfield(model.Facilities, name)
        fac = model.Facilities.(name);

        if isfield(fac, 'Type') && strcmp(fac.Type, 'M')
            [newX, newY] = randomWalk(rect.x, rect.y, ...
                                      model.site_width, model.site_height, ...
                                      rect.width, rect.height);
            sol.Stage(stageIdx).Facilities.(name).x = newX;
            sol.Stage(stageIdx).Facilities.(name).y = newY;
        end

    % Materials
    elseif isfield(model.materials, name)
        [newX, newY] = randomWalk(rect.x, rect.y, ...
                                  model.site_width, model.site_height, ...
                                  rect.width, rect.height);
        sol.Stage(stageIdx).Materials.(name).x = newX;
        sol.Stage(stageIdx).Materials.(name).y = newY;
    end
end
function d = rectDistance(A, B)
    % A و B ساختار با فیلدهای x, y, Lx, Ly هستند
    dx = max(0, max(B.x - (A.x + A.width), A.x - (B.x + B.width)));
    dy = max(0, max(B.y - (A.y + A.height), A.y - (B.y + B.height)));
    d = sqrt(dx^2 + dy^2);
end
