function [x, y] = FindFreePositionGreedy(gridMap, Lx, Ly, blockedZones)
    [H, W] = size(gridMap);
    x = -1;
    y = -1;

        % --- 1. نواحی ممنوعه رو در gridMap علامت‌گذاری کن ---
    for i = 1:numel(blockedZones)
        bx = blockedZones(i).x;
        by = blockedZones(i).y;
        bw = blockedZones(i).width;
        bh = blockedZones(i).height;

        % محاسبه اندیس‌ها با اطمینان از اینکه از مرز بیرون نزنیم
        x1 = max(1, bx);
        y1 = max(1, by);
        x2 = min(W, bx + bw - 1);
        y2 = min(H, by + bh - 1);

        gridMap(y1:y2, x1:x2) = 1;
    end

    % از بالا به پایین و چپ به راست بگرد
    for row = 1:H - Ly + 1
        for col = 1:W - Lx + 1
            sub = gridMap(row:row+Ly-1, col:col+Lx-1);
            if all(sub(:) == 0)
                x = col;
                y = row;
                return;
            end
        end
    end
end
