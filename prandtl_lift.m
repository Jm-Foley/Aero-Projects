clear
clc

filename = input('ascii filename: ','s');
data = load(filename);
disp ('wing geometric parameters: ')
S = input('wing area: ');
lambda = input('taper ratio: ');
AR = input('aspect ratio: ');
ttwist = deg2rad(input('tip twist angle: '));

disp ('environment parameters: ')
h = input('altitude: ');
Vinf = input('freestream velocity: ');
b = sqrt(AR*S);
bhalf = 1/2*b;

Croot = 2*(S/(b*(1+lambda)));

[~, ~, ~, rho] = atmosisa(h);

AoA = deg2rad(data(:,1));
Cl = data(:,2);
Cd = data(:,3);



function w = hshoe(y,y1,Gamma)  
% gamma column to sort dimension issue
    Gamma = Gamma(:);
    w = Gamma .* (1./(y - y1(:,1)) + 1./(y1(:,2) - y)) / (4*pi);  
end

% create function
function [CD,CL,Cdint,Clint,y,c] = wingfunction(AoAc,N,ttwist,lambda,Croot,Vinf,b,Cl,Cd,AoA,S)
    widthI = b / N;
    y = linspace(-b/2 + widthI/2, b/2 - widthI/2, N);
    y1 = [(y - widthI/2)', (y + widthI/2)'];
    c = Croot * (1 - (1 - lambda) * abs(y)/(b/2));  
    eps = ttwist * (abs(y)/(b/2));  

    % initial conditions for iteration
    circ = ((4*Vinf*S*AoAc/b).*sqrt(1-(y/(b/2)).^2));
    circNew = circ+1;
    tol = 1e-6;
    D = 0.01;
    W = zeros(1,N);
    % iteration to find real circulation
    while max(abs(circ-circNew))>tol
        for i = 1:N
            
            W(i) = sum(hshoe(y(i), y1, circ));
            
        end
        % find new effective AoA
        AoAe = eps + atan2((Vinf*sin(AoAc)-W),(Vinf*cos(AoAc)));

        % interpolate at new AoA
        Clint = interp1(AoA,Cl,AoAe,'linear');
        Cdint = interp1(AoA,Cd,AoAe,'linear');
        
        % kutta-joukowski for new circulation
        circNew = (0.5 * Vinf * c .* Clint);

        % apply damping ratio
        circ = circ + D*(circNew - circ);
        
    end
    % find change in AoA
    dAoA = AoAe - eps - AoAc;
    Ycl = c.*Clint.*cos(dAoA);
    Ycd = c.*(Cdint - Clint.*(sin(dAoA)));
    % apply simpsons 1/3
    CL = (1/S) * (widthI/3) * (Ycl(1) + 4*sum(Ycl(2:2:end)) + 2*sum(Ycl(3:2:end-1)) + Ycl(end));
    CD = (1/S) * (widthI/3) * (Ycd(1) + 4*sum(Ycd(2:2:end)) + 2*sum(Ycd(3:2:end-1)) + Ycd(end));
end

mode = input('choose mode; mode 1 (single angle of attack and number of sections), mode 2 (fixed number of panels with range of AoA), mode 3 (fixed AoA with varying number of panels): ');
while true
    if mode == 1
        % user input for mode 1, 1 AoA chosen and 1 number of sections
        AoAc = deg2rad(input('chosen AoA: '));
        N = input ('number of sections: ');
        [CD,CL,Cdint,Clint,y,c] = wingfunction(AoAc,N,ttwist,lambda,Croot,Vinf,b,Cl,Cd,AoA,S);
        % wing tip constraints
        Clint(1) = 0;
        Clint(end) = 0;
        % find lift with cl, S, rho, and Vinf
        Lint = 1/2*rho*(Vinf^2)*c.*Clint;
        % plots for mode 1
        figure()
        plot(y, Lint, '-r', 'LineWidth', 1.5)
        xlabel('Spanwise position y (m)')
        ylabel('Lift per unit span (N/m)')
        title(['Lift distribution, AoA = ',num2str(rad2deg(AoAc)),'deg, N =',num2str(N)])
        % save/ export plot
        saveas(gcf, ['liftperunitspan,AoA' num2str(rad2deg(AoAc)) 'deg,N' num2str(N) '.png'])
        % save text file
        diary(['M1results,AoA' num2str(rad2deg(AoAc)) 'deg,N' num2str(N) '.txt'])
        disp('========================================')
        disp('Wing Analysis Results - Mode 1')
        disp('========================================')
        disp(['AoA = ' num2str(rad2deg(AoAc)) ' deg'])
        disp(['N = ' num2str(N) ' panels'])
        disp(['Altitude = ' num2str(h) ' m'])
        disp(['Vinf = ' num2str(Vinf) ' m/s'])
        disp('----------------------------------------')
        disp(['3D Lift Coefficient CL: ' num2str(CL)])
        disp(['3D Drag Coefficient CD: ' num2str(CD)])
        disp('========================================')
        diary off

        break
    elseif mode == 2
        % user input for mode 2, chosen N and range of AoA
        N = input ('number of sections: ');
        AoAc = deg2rad(linspace(-20,30,20));
        % initialise 3d Cd and Cl, final Cd and Cl in 2d for each section
        CD = zeros(size(AoAc));
        CL = zeros(size(AoAc));
        Cdint = zeros(length(AoAc), N);
        Clint = zeros(length(AoAc), N);
        y = zeros(N);
        colours = lines(length(AoAc));
        figure()
        for i = 1:length(AoAc) % iterate through function with each AoA configuration
            [CD(i),CL(i),Cdint(i,:),Clint(i,:),y,c] = wingfunction(AoAc(i),N,ttwist,lambda,Croot,Vinf,b,Cl,Cd,AoA,S); 
        end
        % 3d cl curve
        figCl = figure();
        AoAdeg = rad2deg(AoAc);
        p = polyfit(AoAdeg, CL, 5);
        AoAsmooth = linspace(min(AoAdeg), max(AoAdeg), 500);
        CLsmooth = polyval(p, AoAsmooth);
        plot(AoAdeg, CL, 'xr', 'LineWidth', 1.5)
        hold on
        plot(AoAsmooth, CLsmooth, '-b', 'LineWidth', 1.5)
        title('3D lift coefficient against AoA')
        xlabel('AoA (degrees)')
        ylabel('3D CL')
        legend('Data points', 'Fitted curve', 'Location', 'eastoutside')
        hold off
        % save figCl
        saveas(figCl, ['CLvsAoA,N' num2str(N) '.png'])
        % 3d cd curve
        figCd = figure();
        p1 = polyfit(AoAdeg, CD, 5);
        CDsmooth = polyval(p1, AoAsmooth);
        plot(AoAdeg, CD, 'xr', 'LineWidth', 1.5)
        hold on
        plot(AoAsmooth, CDsmooth, '-b', 'LineWidth', 1.5)
        title('3D drag coefficient against AoA')
        xlabel('AoA (degrees)')
        ylabel('3D CD')
        legend('Data points', 'Fitted curve', 'Location', 'eastoutside')
        hold off
        % save figCd
        saveas(figCd, ['CDvsAoA,N' num2str(N) '.png'])
        % textfile
        diary(['M2results,N' num2str(N) '.txt'])
        disp('========================================')
        disp('Wing Analysis Results - Mode 2')
        disp('========================================')
        disp(['N = ' num2str(N) ' panels'])
        disp(['Altitude = ' num2str(h) ' m'])
        disp(['Vinf = ' num2str(Vinf) ' m/s'])
        disp('----------------------------------------')
        disp('AoA (deg)    CL        CD')
        disp('----------------------------------------')
        for i = 1:length(AoAc)
            disp([num2str(rad2deg(AoAc(i))) '    ' num2str(CL(i)) '    ' num2str(CD(i))])
        end
        disp('========================================')
        diary off
        break
    elseif mode == 3
        % take input of AoA, and set range of N
        AoAc = deg2rad(input('chosen AoA: '));
        N = [10,20,50,100];
        % initialise 3d Cd and Cl, and cd and cl for each 2d section
        CD = zeros(size(N));
        CL = zeros(size(N));
        Clint = cell(length(N),1);
        Cdint = cell(length(N),1);
        y = cell(length(N),1);
        c = cell(length(N),1);
        colours = lines(length(N));
        figure()
        for i = 1:length(N) % iterate through wingfunction with each number of sections N(i)
            [CD(i),CL(i),Cdint{i},Clint{i},y{i},c{i}] = wingfunction(AoAc,N(i),ttwist,lambda,Croot,Vinf,b,Cl,Cd,AoA,S);
            Clint{i}(1) = 0;
            Clint{i}(end) = 0;
            % find lift with cl, S, rho, and Vinf
            Lint{i} = 1/2*rho*(Vinf^2)*c{i}.*Clint{i};
            plot(y{i},Lint{i},'-','Color', colours(i,:), 'LineWidth', 1.5)
            hold on
            legendEntries{i} = (['AoA = ',num2str(rad2deg(AoAc)),'deg, N =',num2str(N(i))]);
        end
        legend(legendEntries, 'Location', 'eastoutside')
        title(['Lift per unit span, AoA =',num2str(rad2deg(AoAc)),' N: 10-100'])
        xlabel('Spanwise position y (m)')
        ylabel('Lift per unit span (N/m)')
        % save plot
        saveas(gcf, ['liftperunitspan,AoA' num2str(rad2deg(AoAc)) 'deg,rangeN.png'])
        % text file
        diary(['M3results,AoA' num2str(rad2deg(AoAc)) 'deg.txt'])
        disp('========================================')
        disp('Wing Analysis Results - Mode 3')
        disp('========================================')
        disp(['AoA = ' num2str(rad2deg(AoAc)) ' deg'])
        disp(['Altitude = ' num2str(h) ' m'])
        disp(['Vinf = ' num2str(Vinf) ' m/s'])
        disp('----------------------------------------')
        disp('N panels    CL        CD')
        disp('----------------------------------------')
        for i = 1:length(N)
            disp([num2str(N(i)) '    ' num2str(CL(i)) '    ' num2str(CD(i))])
        end
        disp('========================================')
        diary off
        break
    else
        mode = input('please choose valid mode (input 1 2 or 3): ');
    end
end


