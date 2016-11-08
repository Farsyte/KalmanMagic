function fig = plot_zxe(ttl, z, H, x, v, z_, x_, fig)

    e = v;
    e(e<0) = 0;
    e = sqrt(e);

    [Nz Nx] = size(H);

    z_e   = H*(x);

    fig = figure('Position', [1920 22 1280 1002],'Name',ttl);

    spr = Nz+Nx;
    spc = 3;
    spi = 0;

    % Change in plans.
    % Don't plot est vs act,
    % plot the difference.

    z   = z   - z_;
    x   = x   - x_;
    z_e = z_e - z_;

    for j=1:Nz

        tlo = 1; thi = 10; t = tlo:thi;

        spi = spi+1; subplot(spr,spc,spi);
        plot([tlo thi],[0 0],'--k');
        hold on;
        plot(t,z(j,t),'.r');
        plot(t,z_e(j,t),'-b');
        zr = max(abs(z_e(j,t)));
        axis([tlo thi -zr +zr]);

        tlo = 10; thi = 100; t = tlo:thi;

        spi = spi+1; subplot(spr,spc,spi);
        plot([tlo thi],[0 0],'--k');
        hold on;
        plot(t,z(j,t),'.r');
        plot(t,z_e(j,t),'-b');
        zr = max(abs(z_e(j,t)));
        axis([tlo thi -zr +zr]);

        tlo = 100; thi = 1000; t = tlo:thi;

        spi = spi+1; subplot(spr,spc,spi);
        plot([tlo thi],[0 0],'--k');
        hold on;
        plot(t,z(j,t),'.r');
        plot(t,z_e(j,t),'-b');
        zr = max(abs(z_e(j,t)));
        axis([tlo thi -zr +zr]);

    end

    for i=1:Nx

        tlo = 1; thi = 10; t = tlo:thi;
        spi = spi+1; subplot(spr,spc,spi);
        plot([tlo thi],[0 0],'--k');
        hold on;
        plot(t,x(i,t)+e(i,t),'-m');
        plot(t,x(i,t)-e(i,t),'-m');
        plot(t,x(i,t)+3*e(i,t),'--c');
        plot(t,x(i,t)-3*e(i,t),'--c');
        zr = max(abs(x(i,t))) + max(abs(e(i,t)));
        axis([tlo thi -zr +zr]);

        tlo = 10; thi = 100; t = tlo:thi;
        spi = spi+1; subplot(spr,spc,spi);
        plot([tlo thi],[0 0],'--k');
        hold on;
        plot(t,x(i,t)+e(i,t),'-m');
        plot(t,x(i,t)-e(i,t),'-m');
        plot(t,x(i,t)+3*e(i,t),'--c');
        plot(t,x(i,t)-3*e(i,t),'--c');
        zr = max(abs(x(i,t))) + max(abs(e(i,t)));
        axis([tlo thi -zr +zr]);

        tlo = 100; thi = 1000; t = tlo:thi;
        spi = spi+1; subplot(spr,spc,spi);
        plot([tlo thi],[0 0],'--k');
        hold on;
        plot(t,x(i,t)+e(i,t),'-m');
        plot(t,x(i,t)-e(i,t),'-m');
        plot(t,x(i,t)+3*e(i,t),'--c');
        plot(t,x(i,t)-3*e(i,t),'--c');
        zr = max(abs(x(i,t))) + max(abs(e(i,t)));
        axis([tlo thi -zr +zr]);

    end

    % for i=1:Nx
    %
    %   spi = spi+1; subplot(spr,spc,spi);
    %   plot(t,e(i,t),'-m');
    %   axis([1 10 0 max(e(i,t))*1.1]);
    %
    %   spi = spi+1; subplot(spr,spc,spi);
    %   plot(10:100,e(i,10:100),'-m');
    %   axis([10 100 0 max(e(i,10:100))*1.1]);
    %
    %   spi = spi+1; subplot(spr,spc,spi);
    %   plot(100:1000,e(i,100:1000),'-m');
    %   axis([100 1000 0 max(e(i,100:1000))*1.1]);
    %
    % end
