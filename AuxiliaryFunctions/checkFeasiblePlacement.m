function ok = checkFeasiblePlacement(newRect, existingRects, constraints)
    ok = true;
    for k = 1:numel(existingRects)
        A = newRect;
        B = existingRects(k);

        % فاصله
        d = rectDistance(A, B);  

        % Safety Min constraints
        for s = 1:size(constraints.safetyMinPairs,1)
            if (strcmp(B.name,constraints.safetyMinPairs{s,1}) && strcmp(A.name,constraints.safetyMinPairs{s,2})) || ...
               (strcmp(B.name,constraints.safetyMinPairs{s,2}) && strcmp(A.name,constraints.safetyMinPairs{s,1}))
                if d < constraints.safetyMinPairs{s,3}
                    ok = false; return;
                end
            end
        end
        
        % Operational Max constraints (مرکز به مرکز)
        for s = 1:size(constraints.operationalMaxPairs,1)
            if (strcmp(B.name,constraints.operationalMaxPairs{s,1}) && strcmp(A.name,constraints.operationalMaxPairs{s,2})) || ...
               (strcmp(B.name,constraints.operationalMaxPairs{s,2}) && strcmp(A.name,constraints.operationalMaxPairs{s,1}))
                % فاصله مرکز به مرکز
                centerA = [A.x + A.width/2, A.y + A.height/2];
                centerB = [B.x + B.width/2, B.y + B.height/2];
                d_center = norm(centerA - centerB);
                if d_center > constraints.operationalMaxPairs{s,3}
                    ok = false; return;
                end
            end
        end

        % Operational Min constraints
        for s = 1:size(constraints.operationalMinPairs,1)
            if (strcmp(B.name,constraints.operationalMinPairs{s,1}) && strcmp(A.name,constraints.operationalMinPairs{s,2})) || ...
               (strcmp(B.name,constraints.operationalMinPairs{s,2}) && strcmp(A.name,constraints.operationalMinPairs{s,1}))
                if d < constraints.operationalMinPairs{s,3}
                    ok = false; return;
                end
            end
        end
    end
end
function d = rectDistance(A, B)
    % A و B باید struct با فیلدهای x, y, width, height باشن
    dx = max(0, max(B.x - (A.x + A.width), A.x - (B.x + B.width)));
    dy = max(0, max(B.y - (A.y + A.height), A.y - (B.y + B.height)));
    d_edge = sqrt(dx^2 + dy^2);

    % فاصله مرکز به مرکز هم
    centerA = [A.x + A.width/2, A.y + A.height/2];
    centerB = [B.x + B.width/2, B.y + B.height/2];
    d_center = norm(centerA - centerB);

    % برگردوندن هر دو فاصله (می‌تونی انتخاب کنی کدوم استفاده بشه)
    d = min(d_edge, d_center);   % یا فقط d_center اگه بخوای
end

