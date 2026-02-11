function sol = ConstraintRepair2(sol, model, it)

    % برای جلوگیری از تغییرات شدید در نسل‌های اولیه، می‌توان ترمیم را به نسل‌های پایانی محدود کرد
    if it < model.MaxIt / 2 
        return;
    end

    MAX_REPAIR_ITERATIONS = 5; % برای جلوگیری از حلقه بی‌نهایت

    for iter = 1:MAX_REPAIR_ITERATIONS
        isRepaired = false; % آیا در این دور تکرار، ترمیمی انجام شد؟

        for i = 1:numel(sol.Stage)
            
            % --- Safety (Min Distance Violations) ---
            for s = 1:size(model.Constraints.safetyMinPairs, 1)
                A_name = model.Constraints.safetyMinPairs{s, 1};
                B_name = model.Constraints.safetyMinPairs{s, 2};
                d_min = model.Constraints.safetyMinPairs{s, 3};
                
                [rectA, rectB] = getRects(sol.Stage(i).occupiedRects, A_name, B_name);
                if isempty(rectA) || isempty(rectB), continue; end

                dist = rectDistance(rectA, rectB);

                if dist < d_min
                    deficit = d_min - dist; % مقدار کمبود فاصله
                    [sol, repaired_now] = intelligentMove(sol, model, i, A_name, B_name, 'separate', deficit);
                    isRepaired = isRepaired || repaired_now;
                end
            end

            % --- Operational (Max Distance Violations) ---
            for s = 1:size(model.Constraints.operationalMaxPairs, 1)
                A_name = model.Constraints.operationalMaxPairs{s, 1};
                B_name = model.Constraints.operationalMaxPairs{s, 2};
                d_max = model.Constraints.operationalMaxPairs{s, 3};

                [rectA, rectB] = getRects(sol.Stage(i).occupiedRects, A_name, B_name);
                if isempty(rectA) || isempty(rectB), continue; end
                
                dist = rectDistance(rectA, rectB);

                if dist > d_max
                    surplus = dist - d_max; % مقدار فاصله اضافی
                    [sol, repaired_now] = intelligentMove(sol, model, i, A_name, B_name, 'converge', surplus);
                    isRepaired = isRepaired || repaired_now;
                end
            end
        end

        % اگر در یک دور کامل هیچ ترمیمی انجام نشد، یعنی به پایداری رسیده و می‌توان خارج شد
        if ~isRepaired
            break;
        end
    end
end