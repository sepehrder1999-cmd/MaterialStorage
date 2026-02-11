function feasible = checkFeasibility(sol, model, stageIdx)
    feasible = true;
    stage = sol.Stage(stageIdx);

    % دسترسی سریع به همه فسیلیتی‌ها و مصالح این مرحله
    items = struct();
    if isfield(stage, 'Facilities') && isstruct(stage.Facilities)
        facilities = fieldnames(stage.Facilities);
    else
        facilities = {};
    end
        for i = 1:length(facilities)
            name = facilities{i};
            item = stage.Facilities.(name);
            items.(name) = item;
        end
    if isfield(stage, 'Materials') && isstruct(stage.Materials)
        materials = fieldnames(stage.Materials);
    else
        materials = {};
    end
        for i = 1:length(materials)
            name = materials{i};
            item = stage.Materials.(name);
            items.(name) = item;
        end

    % Helper برای گرفتن مرکز
    getCenter = @(r) [r.x + r.Lx/2, r.y + r.Ly/2];

    % بررسی safetyMinPairs
    pairs = model.Constraints.safetyMinPairs;
    for i = 1:size(pairs,1)
        a = pairs{i,1};
        b = pairs{i,2};
        minDist = pairs{i,3};

        if isfield(items, a) && isfield(items, b)
            ra = items.(a);
            rb = items.(b);
            ca = getCenter(ra);
            cb = getCenter(rb);
            d = norm(ca - cb);

            if d < minDist
                feasible = false;
                return;
            end
        end
    end

    % بررسی operationalMaxPairs
    pairs = model.Constraints.operationalMaxPairs;
    for i = 1:size(pairs,1)
        a = pairs{i,1};
        b = pairs{i,2};
        maxDist = pairs{i,3};

        if isfield(items, a) && isfield(items, b)
            ra = items.(a);
            rb = items.(b);
            ca = getCenter(ra);
            cb = getCenter(rb);
            d = norm(ca - cb);

            if d > maxDist
                feasible = false;
                return;
            end
        end
    end

    % بررسی operationalMinPairs
    pairs = model.Constraints.operationalMinPairs;
    for i = 1:size(pairs,1)
        a = pairs{i,1};
        b = pairs{i,2};
        minDist = pairs{i,3};

        if isfield(items, a) && isfield(items, b)
            ra = items.(a);
            rb = items.(b);
            ca = getCenter(ra);
            cb = getCenter(rb);
            d = norm(ca - cb);

            if d < minDist
                feasible = false;
                return;
            end
        end
    end
end
