function plotTest(concentrationsM0, concentrationsM1, concentrationsM2, gene)
% This function will generate a plot with expected and our sfc curves of all three morphogen conditions

datag1 = concentrationsM0(1:3:96);
datag2 = concentrationsM0(2:3:96);
datag3 = concentrationsM0(3:3:96);

 datag1M1 = concentrationsM1(1:3:96);
 datag2M1 = concentrationsM1(2:3:96);
 datag3M1 = concentrationsM1(3:3:96);

datag1M2 = concentrationsM2(1:3:96);
datag2M2 = concentrationsM2(2:3:96);
datag3M2 = concentrationsM2(3:3:96);

if (min(datag1) ~= max(datag1))
    datag1 = maptorange(datag1, [min(datag1), max(datag1)], [0, 2]);
end
if (min(datag2) ~= max(datag2))
    datag2 = maptorange(datag2, [min(datag2), max(datag2)], [0, 2]);
end

if (min(datag3) ~= max(datag3))
    datag3 = maptorange(datag3, [min(datag3), max(datag3)], [0, 2]);
end
if (min(datag1M1) ~= max(datag1M1))
    datag1M1 = maptorange(datag1M1, [min(datag1M1), max(datag1M1)], [0, 2]);
end
if (min(datag2M1) ~= max(datag2M1))
    datag2M1 = maptorange(datag2M1, [min(datag2M1), max(datag2M1)], [0, 2]);
end

if (min(datag3M1) ~= max(datag3M1))
    datag3M1 = maptorange(datag3M1, [min(datag3M1), max(datag3M1)], [0, 2]);
end
if (min(datag1M2) ~= max(datag1M2))
    datag1M2 = maptorange(datag1M2, [min(datag1M2), max(datag1M2)], [0, 2]);
end
if (min(datag2M2) ~= max(datag2M2))
    datag2M2 = maptorange(datag2M2, [min(datag2M2), max(datag2M2)], [0, 2]);
end

if (min(datag3M2) ~= max(datag3M2))
    datag3M2 = maptorange(datag3M2, [min(datag3M2), max(datag3M2)], [0, 2]);
end

sfcM2 = [0.0221047858304046,0.0352278449085936,0.0554442306835930,0.0856265102987635,0.128650047912544,0.186128995225380,0.256703382176614,0.335146725042178,0.413861738245364,0.486504776127235,0.551283143017491,0.612319565692932,0.680641940101599,0.779779499174902,0.973241271495922,1.55605548046255,28.9050764890140,-1.05925239429263,-0.395170903303485,-0.195716342758060,-0.107459416457579,-0.0619694996030804,-0.0366855347004706,-0.0220424831085947,-0.0133598512571004,-0.00813946086272128,-0.00497449900841864,-0.00304599569637634,-0.00186729685926280,-0.00114552835080721,-0.000703051693030268];
sfc = [0.00761819493947260,0.0102073826846778,0.0136422707833118,0.0181726316482285,0.0241022131401909,0.0317858937828222,0.0416147301162031,0.0539824030656485,0.0692276737132812,0.0875530646830387,0.108931735080854,0.133030383231278,0.159187354505148,0.186478859382782,0.213875807016407,0.240451958286658,0.265579469893321,0.289060796822430,0.311188626135617,0.332771962823697,0.355200928276759,0.380659435820751,0.412690396614078,0.457638499437845,0.528692071310469,0.659855717272589,0.977415600284468,2.70501870948564];
x2 = [0.973910000000000,1.87830000000000,3.47830000000000,5.21740000000000,7.72170000000000,10.2260000000000,11.0610000000000,13.1480000000000,13.8430000000000,15.6520000000000,18.4350000000000,21.5650000000000,23.0260000000000,23.8610000000000,24.9040000000000,25.3220000000000,25.8780000000000,26.2260000000000,26.5040000000000,26.7130000000000,26.8520000000000];
y2 = [0.00251890000000000,0.00251890000000000,0.00503780000000000,0.0125940000000000,0.0428210000000000,0.0680100000000000,0.0856420000000000,0.130980000000000,0.138540000000000,0.161210000000000,0.216620000000000,0.282120000000000,0.317380000000000,0.365240000000000,0.418140000000000,0.486150000000000,0.604530000000000,0.717880000000000,0.851390000000000,0.947100000000000,1.00500000000000];
for i = 1%:5%:100
    x=1:32;
    y=1:28;
    y2 = 1:17;
    figure(i)
if strcmp(gene, 'g1')==1
    plot(x,datag1,'g','DisplayName','g1'); % plot of g3
    hold on
     plot(x,datag1M1,'g*','DisplayName','g1M1'); % plot of g3
     hold on
     plot(x,datag1M2,'-','DisplayName','g1M2'); % plot of g3
     hold on
end
if strcmp(gene, 'g2')==1
    plot(x,datag2,'g','DisplayName','g2'); % plot of g3
    hold on
     plot(x,datag2M1,'g*','DisplayName','g2M1'); % plot of g3
     hold on
     plot(x,datag2M2,'-','DisplayName','g2M2'); % plot of g3
     hold on
end
if strcmp(gene, 'g3')==1
    plot(x,datag3,'g','DisplayName','g3'); % plot of g3
    hold on
     plot(x,datag3M1,'g*','DisplayName','g3M1'); % plot of g3
     hold on
     plot(x,datag3M2,'-','DisplayName','g3M2'); % plot of g3
     hold on
end
    plot(y,sfc(1:28),'m','DisplayName','sfc')
    hold on
    plot(y2,sfcM2(1:17),'m','DisplayName','sfcM2')
    axis([1 40 0 1]);
    legend
    %hold on
  %   plot(x2,y2,'b');
  %  plot(x,z,'*r','DisplayName','M'); % plot of Morphogen
end
end