function PlotSolution(sol, model)

    sol = ParseSolution(sol, model);
    figure(1);          % همیشه از شکل شماره 1 استفاده کن
    clf;                % پاک‌سازی شکل برای رسم جدید
    T = numel(sol.Stage);
    
    for i = 1:T
        subplot(1, T, i);
        hold on;
        axis equal;
        title(['Stage ' num2str(i) ' Layout']);
        xlabel('Width (m)');
        ylabel('Height (m)');
        xlim([0, model.site_width]);
        ylim([0, model.site_height]);

        % مرز سایت
        rectangle('Position', [0 0 model.site_width model.site_height], ...
          'EdgeColor', 'k', 'LineWidth', 1);

    
        % Buildings
        if isfield(sol.Stage(i), 'Buildings')
            all_bnames = fieldnames(sol.Stage(i).Buildings);
            bnames = all_bnames(startsWith(all_bnames, 'B'));  % فقط ساختمان‌های B            
            for j = 1:numel(bnames)
                bname = bnames{j};
                b = sol.Stage(i).Buildings.(bname);
                delta = 0.75;
                rectangle('Position', [b.x+delta, b.y+delta, b.Lx-2*delta, b.Ly-2*delta], ...
                    'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k', 'LineWidth', 1.5);
                text(b.x + b.Lx/2, b.y + b.Ly/2, bname, ...
                    'HorizontalAlignment', 'center', 'FontSize', 10);
            end
        end
    
        % Gate
        if isfield(model.Facilities, 'G')
            g = model.Facilities.G;
            rectangle('Position', [g.x, g.y, g.Lx, g.Ly], ...
                'EdgeColor', 'k', 'LineWidth', 1);
        end
    
        % Facilities
        if isfield(sol.Stage(i), 'Facilities')
            fnames = fieldnames(sol.Stage(i).Facilities);
            for j = 1:numel(fnames)
                fname = fnames{j};
                f = sol.Stage(i).Facilities.(fname);
                rectangle('Position', [f.x+delta, f.y+delta, f.Lx-2*delta, f.Ly-2*delta], ...
                    'EdgeColor', 'g', 'FaceColor', [0.8 1 0.8], 'LineWidth', 1.5);
                text(f.x + f.Lx/2, f.y + f.Ly/2, fname, ...
                    'HorizontalAlignment', 'center', 'FontSize', 8);
            end
        end
    
        % Materials
        if isfield(sol.Stage(i), 'Materials')
            mnames = fieldnames(sol.Stage(i).Materials);
            for j = 1:numel(mnames)
                mid = mnames{j};
                m = sol.Stage(i).Materials.(mid);
                rectangle('Position', [m.x+delta, m.y+delta, m.Lx-2*delta, m.Ly-2*delta], ...
                    'EdgeColor', 'r', 'FaceColor', [1 0.8 0.8], 'LineWidth', 1.5);
                text(m.x + m.Lx/2, m.y + m.Ly/2, mid, ...
                    'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'r');
            end
        end
        % ===== نمایش ناحیه exclusion =====
        if isfield(model, 'Constraints') && isfield(model.Constraints, 'exclusion')
            exc = model.Constraints.exclusion;
            if ismember(i, exc.activeStages)
                x = exc.x;
                y = exc.y;
                w = exc.Lx;
                h = exc.Ly;
        
                % رسم ناحیه با خط‌چین مشکی
                rectangle('Position', [x y w h], 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
        
                % برچسب exc.
                text(x + w/2, y + h/2, 'exc.', 'FontSize', 12, 'Color', 'k', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            end
        end

    end
end
