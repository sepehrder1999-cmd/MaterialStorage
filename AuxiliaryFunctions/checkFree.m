function free = checkFree(gridMap, x, y, w, h, blockedZones)

    H = size(gridMap,1);
    W = size(gridMap,2);
    
    if x < 1 || y < 1 || x+w-1 > W -0.5|| y+h-1 > H -0.5
        free = false;
        return;
    end
    sub = gridMap(y:y+h-1, x:x+w-1);
    free = all(sub(:) == 0);

    if ~free  || isempty(blockedZones)
        return;
    end

    for i = 1:numel(blockedZones)
        bx = blockedZones(i).x;
        by = blockedZones(i).y;
        bw = blockedZones(i).width;
        bh = blockedZones(i).height;
    
        if rectsOverlap([x, y, w, h], [bx, by, bw, bh])
            free = false;
            return;
        end
    end
end

function overlap = rectsOverlap(r1, r2)
    x1 = r1(1); y1 = r1(2); w1 = r1(3); h1 = r1(4);
    x2 = r2(1); y2 = r2(2); w2 = r2(3); h2 = r2(4);
    
    overlap = ~(x1+w1 < x2 || x2+w2 < x1 || y1+h1 < y2 || y2+h2 < y1);
end
