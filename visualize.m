% VISUALIZING THE GRAPHS
parameters = newParameters2;
digraphsFinal = randDigraphs2000;
params = 3;
digraphs = 1;
%sfc ofr M0 and M1
sfc1 = [0.00760273448820800,0.0101864566941910,0.0136139222327332,0.0181341842428832,0.0240499942911999,0.0317148437179244,0.0415178479478744,0.0538499604451974,0.0690461005234863,0.0873033787501381,0.108587347833183,0.132554059276025,0.158527016470830,0.185561770269678,0.212600363327769,0.238675783939125,0.263101534166405,0.285593520864504,0.306312555044299,0.325858585715879,0.345273694231898,0.366125561978866,0.390776157991274,0.423056762731367,0.469965559726712,0.546503189197296,0.693131717513635,1.07359413237585,4.00601130890706,-1.50471828461896,-0.527969279625462];
%sfc for M2
sfc2 = [0.0221047858304046,0.0352278449085936,0.0554442306835930,0.0856265102987635,0.128650047912544,0.186128995225380,0.256703382176614,0.335146725042178,0.413861738245364,0.486504776127235,0.551283143017491,0.612319565692932,0.680641940101599,0.779779499174902,0.973241271495922,1.55605548046255,28.9050764890140,-1.05925239429263,-0.395170903303485,-0.195716342758060,-0.107459416457579,-0.0619694996030804,-0.0366855347004706,-0.0220424831085947,-0.0133598512571004,-0.00813946086272128,-0.00497449900841864,-0.00304599569637634,-0.00186729685926280,-0.00114552835080721,-0.000703051693030268];

% mode = 1 means it runs on Michaelis-Menten Function, 2 means run on Sigmoid function
mode = 1;
% parameters for the MM function or Sigmoid function
param1 = 5;
param2 = 5;

% 12 parameters and digraphs used
para = 1;
digraph =1;




        genes_M0 = modelEuler_v3(zeros(1,96),parameters(para,:),digraphsFinal{digraph}, mode, param1,param2);    
        mapped_M0 = map_M_torange(genes_M0);
        % visualizing M0
        g1_M0 = mapped_M0(1,:);
        g2_M0 = mapped_M0(2,:);
        g3_M0 = mapped_M0(3,:);
        length = size(g1_M0);
        lengthsfc = size(sfc1); 
        figure(1);
        title("M0");
        xlabel("cell");
        ylabel("sfc");
        x = 1:length;
        x_sfc = 1:lengthsfc;
        plot(x,g1_M0);
        hold on;
        plot(x,g2_M0);
        plot(x,g3_M0);
        plot(x_sfc, sfc1);
        
        
        
        
        
        
      
        
        genes_M1 = modelEulerM1_v3(zeros(1,96),parameters(para,:),digraphsFinal{digraph}, mode, param1,param2);
        mapped_M1 = map_M_torange(genes_M1);
        % visualize M1
        g1_M1 = mapped_M1(1,:);
        g2_M1 = mapped_M1(2,:);
        g3_M1 = mapped_M1(3,:);
        length = size(g1_M1);
        figure(2);
        title("M1");
        xlabel("cell");
        ylabel("sfc");
        x_sfc = 1:lengthsfc;
        plot(x,g1_M1);
        hold on;
        plot(x,g2_M1);
        plot(x,g3_M1);
        plot(x_sfc, sfc1);
        
        
        
        genes_M2 = modelEulerM2_v3(zeros(1,96),parameters(para,:),digraphsFinal{digraph}, mode, param1,param2);
        mapped_M2 = map_M2_torange(genes_M2);
        % visualize M1
        g1_M2 = mapped_M2(1,:);
        g2_M2 = mapped_M2(2,:);
        g3_M2 = mapped_M2(3,:);
        length = size(g1_M2);
        lengthsfc = size(sfc2); 
        figure(3);
        title("M2");
        xlabel("cell");
        ylabel("sfc");
        x_sfc = 1:lengthsfc;
        plot(x,g1_M2);
        hold on;
        plot(x,g2_M2);
        plot(x,g3_M2);
        plot(x_sfc, sfc2);
        
        
        
%         final_score_G1= (score_M0(1)+score_M1(1)+score_M2(1));
%         final_score_G2= (score_M0(2)+score_M1(2)+score_M2(2));
%         final_score_G3= (score_M0(3)+score_M1(3)+score_M2(3));
    
