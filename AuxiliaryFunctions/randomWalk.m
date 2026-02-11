function [xnew, ynew] = randomWalk(x, y, W, H, w, h)
  dirs = [0 1; 0 -1; 1 0; -1 0]; % بالا، پایین، راست، چپ
    stepSize = max(1, floor(min(W,H)/20)); % اندازه گام قابل تنظیم
    step = dirs(randi(4), :) * stepSize;

    % حرکت و محدود کردن به داخل سایت
    xnew = min(max(x + step(1), 1), W - w + 1);
    ynew = min(max(y + step(2), 1), H - h + 1);
end