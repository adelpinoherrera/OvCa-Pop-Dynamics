function plot_treatment_schedules()
    % Time vector (days, 0–30, 24h steps)
    t = 0:1:30;
    nT = numel(t);

    % Define MTD (arbitrary units; label as MTD in figure)
    MTD = 1.0;

    % Preallocate
    dose_noTx   = zeros(size(t));
    dose_MTD    = MTD * ones(size(t));

    % For demonstration, define a synthetic population trajectory
    % and a threshold to drive AT decisions.
    % Here we just make something plausible; replace with your real model output.
    N = 1e4 * exp(0.15*t) ./ (1 + exp(0.15*t - 3)); % smooth S-shaped growth
    N = N + 0.05*max(N)*sin(0.8*t); % add some wiggles
    N(N<0) = 1e3;

    % Relative change over previous 24h (ΔN/N_prev)
    dN = [0, diff(N)./N(1:end-1)];
    dN(dN<0) = 0; % for simplicity, treat decreases as 0 in this example

    threshold = 0.5; % example threshold (50% change)

    % Adaptive Therapy 1: symmetric dose modulation
    dose_AT1 = zeros(size(t));
    dose_AT1(1) = 0.5*MTD; % start at mid-dose
    for i = 2:nT
        if dN(i) >= threshold
            dose_AT1(i) = min(MTD, dose_AT1(i-1) + 0.25*MTD);
        else
            dose_AT1(i) = max(0, dose_AT1(i-1) - 0.25*MTD);
        end
    end

    % Adaptive Therapy 2: asymmetric with hold
    dose_AT2 = zeros(size(t));
    dose_AT2(1) = 0.5*MTD;
    for i = 2:nT
        if dN(i) >= threshold
            dose_AT2(i) = min(MTD, dose_AT2(i-1) + 0.25*MTD);
        elseif dN(i) <= 0  % population decreased
            dose_AT2(i) = max(0, dose_AT2(i-1) - 0.25*MTD);
        else
            dose_AT2(i) = dose_AT2(i-1); % hold
        end
    end

    % Adaptive Therapy 3: on/off (dose skipping)
    dose_AT3 = zeros(size(t));
    dose_AT3(1) = MTD;
    for i = 2:nT
        if dN(i) >= threshold
            dose_AT3(i) = MTD;
        else
            dose_AT3(i) = 0;
        end
    end

    % Bang-bang: alternate every 24h
    dose_BB = zeros(size(t));
    for i = 1:nT
        if mod(i-1,2) == 0
            dose_BB(i) = MTD;
        else
            dose_BB(i) = 0;
        end
    end

    % Plot
    figure('Color','w','Position',[100 100 900 700]);

    yticks_vals = [0 0.5 1];
    yticks_labels = {'0','0.5×MTD','MTD'};

    % No treatment
    subplot(6,1,1);
    stairs(t, dose_noTx, 'Color',[0.6 0.6 0.6],'LineWidth',1.5);
    ylim([-0.1 1.1]);
    set(gca,'YTick',yticks_vals,'YTickLabel',yticks_labels);
    ylabel('Dose');
    title('No treatment','FontWeight','normal');
    set(gca,'XTick',[]); grid on;

    % Continuous MTD
    subplot(6,1,2);
    stairs(t, dose_MTD, 'Color',[0.85 0.2 0.2],'LineWidth',1.5);
    ylim([-0.1 1.1]);
    set(gca,'YTick',yticks_vals,'YTickLabel',yticks_labels);
    ylabel('Dose');
    title('Continuous MTD','FontWeight','normal');
    set(gca,'XTick',[]); grid on;

    % AT1
    subplot(6,1,3);
    stairs(t, dose_AT1, 'Color',[0.2 0.5 0.8],'LineWidth',1.5);
    ylim([-0.1 1.1]);
    set(gca,'YTick',yticks_vals,'YTickLabel',yticks_labels);
    ylabel('Dose');
    title('Adaptive Therapy 1 (AT1): dose modulation, symmetric','FontWeight','normal');
    set(gca,'XTick',[]); grid on;
    annotation('textbox',[0.62 0.68 0.35 0.08],...
        'String',{'ΔN ≥ threshold → +0.25 MTD', 'ΔN < threshold/↓ → −0.25 MTD'},...
        'FontSize',7,'FitBoxToText','on','EdgeColor','none');

    % AT2
    subplot(6,1,4);
    stairs(t, dose_AT2, 'Color',[0.2 0.7 0.4],'LineWidth',1.5);
    ylim([-0.1 1.1]);
    set(gca,'YTick',yticks_vals,'YTickLabel',yticks_labels);
    ylabel('Dose');
    title('Adaptive Therapy 2 (AT2): dose modulation, asymmetric','FontWeight','normal');
    set(gca,'XTick',[]); grid on;
    annotation('textbox',[0.62 0.53 0.35 0.1],...
        'String',{'ΔN ≥ threshold → +0.25 MTD', 'ΔN < 0 → −0.25 MTD', '0 ≤ ΔN < threshold → hold'},...
        'FontSize',7,'FitBoxToText','on','EdgeColor','none');

    % AT3
    subplot(6,1,5);
    stairs(t, dose_AT3, 'Color',[0.8 0.3 0.85],'LineWidth',1.5);
    ylim([-0.1 1.1]);
    set(gca,'YTick',yticks_vals,'YTickLabel',yticks_labels);
    ylabel('Dose');
    title('Adaptive Therapy 3 (AT3): dose skipping (0 vs MTD)','FontWeight','normal');
    set(gca,'XTick',[]); grid on;
    annotation('textbox',[0.62 0.38 0.35 0.08],...
        'String',{'ΔN ≥ threshold → MTD', 'ΔN < threshold/↓ → 0'},...
        'FontSize',7,'FitBoxToText','on','EdgeColor','none');

    % Bang-bang
    subplot(6,1,6);
    stairs(t, dose_BB, 'Color',[0.9 0.6 0.1],'LineWidth',1.5);
    ylim([-0.1 1.1]);
    set(gca,'YTick',yticks_vals,'YTickLabel',yticks_labels);
    ylabel('Dose');
    xlabel('Time (days)');
    title('Bang-bang skipping (alternate 0/MTD every 24 h)','FontWeight','normal');
    grid on;

    % Overall title
    sgtitle('Treatment Schedules: No Tx, Continuous MTD, AT1–AT3, and Bang-Bang');

    % Optional: save figure
    % saveas(gcf,'treatment_schedules.png');
end