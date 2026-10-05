%% 20250407 Qingyu

function LDEOutLIBy = I_LGNc1_3DSigmoid(L4EUse,L4IUse,L6LibInd)  %#codegen
     
 
    % first we predict which segment(s) we should use by a rough estimation

    a1 = 0.8775840997695923;
    a2 = -1.5060482025146484;
    a3 = 1.9315876960754395;
    L = 368.0876770019531;
    k = 2.8388843536376953;
    x0 = 4.796327114105225;
            
    c1a = -0.3691696524620056;
    c1b = 0.7570117115974426;
    c1c = -17.637189865112305;
    c2a = 0.6164760589599609;
    c2b = 2.69746994972229;
    c2c = 0.46352824568748474;
    c3a = -0.6773039102554321;
    c3b = 1.1515681743621826;
    c3c = -17.21953010559082;
    L1 = 0.2565772533416748;
    k1 = 14.096871376037598;
    x_01 = -0.1621517390012741;
    L2 = 0.2662561535835266;
    k2 = 13.899486541748047;
    x_02 = -1.6365851163864136;
    L3 = 0.27763885259628296;
    k3 = 13.592023849487305;
    x_03 = -0.0037616153713315725;

    b1 = 0.0001*(c1a*L4EUse+c1b*L4IUse+c1c*L6LibInd);
    b2 = 0.0001*(c2a*L4EUse+c2b*L4IUse+c2c*L6LibInd);
    b3 = 0.0001*(c3a*L4EUse+c3b*L4IUse+c3c*L6LibInd);

    b1 = L1./(1+exp(-k1*(b1-x_01)));
    b2 = L2./(1+exp(-k2*(b2-x_02)));
    b3 = L3./(1+exp(-k3*(b3-x_03)));

    h = 0.001*a1*(1+b1).*L4EUse+0.001*a2*(1+b2).*L4IUse +0.1*a3*(1+b3).*L6LibInd;
    LDEOutLIBy_rough =  L./(1+exp(-k*0.1*(h-x0)));

    LDEOutLIBy = zeros(size(LDEOutLIBy_rough));

    % Define masks for different conditions
    idx1 = LDEOutLIBy_rough < 10;
    idx9 = LDEOutLIBy_rough > 20 & LDEOutLIBy_rough < 40;
    idx2 = LDEOutLIBy_rough > 50 & LDEOutLIBy_rough < 90;
    idx3 = LDEOutLIBy_rough > 100 & LDEOutLIBy_rough < 150;
    idx4 = LDEOutLIBy_rough > 160;

    idx5 = LDEOutLIBy_rough >= 10 & LDEOutLIBy_rough <= 20;
    idx6 = LDEOutLIBy_rough >= 40 & LDEOutLIBy_rough <= 50;
    idx7 = LDEOutLIBy_rough >= 90 & LDEOutLIBy_rough <= 100;
    idx8 = LDEOutLIBy_rough >= 150 & LDEOutLIBy_rough <= 160;


    % Direct assignments for known intervals
    LDEOutLIBy(idx1) = I_LGNc1_2_20(L4EUse(idx1), L4IUse(idx1), L6LibInd(idx1));
    LDEOutLIBy(idx2) = I_LGNc1_40_100(L4EUse(idx2), L4IUse(idx2), L6LibInd(idx2));
    LDEOutLIBy(idx3) = I_LGNc1_90_160(L4EUse(idx3), L4IUse(idx3), L6LibInd(idx3));
    LDEOutLIBy(idx4) = I_LGNc1_150_350(L4EUse(idx4), L4IUse(idx4), L6LibInd(idx4));
    LDEOutLIBy(idx9) = I_LGNc1_10_50(L4EUse(idx9), L4IUse(idx9), L6LibInd(idx9));

     % Interpolated values
    LDEOutLIBy(idx5) = ((LDEOutLIBy_rough(idx5)-10)/10) .* I_LGNc1_2_20(L4EUse(idx5), L4IUse(idx5), L6LibInd(idx5)) + ...
                             ((20-LDEOutLIBy_rough(idx5))/10) .* I_LGNc1_10_50(L4EUse(idx5), L4IUse(idx5), L6LibInd(idx5));

    LDEOutLIBy(idx6) = ((LDEOutLIBy_rough(idx6)-40)/10) .* I_LGNc1_10_50(L4EUse(idx6), L4IUse(idx6), L6LibInd(idx6)) + ...
                             ((50-LDEOutLIBy_rough(idx6))/10) .* I_LGNc1_40_100(L4EUse(idx6), L4IUse(idx6), L6LibInd(idx6));

    LDEOutLIBy(idx7) = ((LDEOutLIBy_rough(idx7)-90)/10) .* I_LGNc1_40_100(L4EUse(idx7), L4IUse(idx7), L6LibInd(idx7)) + ...
                             ((100-LDEOutLIBy_rough(idx7))/10) .* I_LGNc1_90_160(L4EUse(idx7), L4IUse(idx7), L6LibInd(idx7));

    LDEOutLIBy(idx8) = ((LDEOutLIBy_rough(idx8)-150)/10) .* I_LGNc1_90_160(L4EUse(idx8), L4IUse(idx8), L6LibInd(idx8)) + ...
                             ((160-LDEOutLIBy_rough(idx8))/10) .* I_LGNc1_150_350(L4EUse(idx8), L4IUse(idx8), L6LibInd(idx8));

end

    function LDEOutLIBy_2_20 = I_LGNc1_2_20(L4ESeg,L4ISeg,L6Ind)
        % parameter for 0-20Hz
        a1 = 1.291521430015564;
        a2 = -2.0586862564086914;
        a3 = 2.970939874649048;
        c1a = 3.3333680629730225;
        c1b = -0.7392067313194275;
        c1c = 0.041597627103328705;
        c2a = 0.6164760589599609;
        c2b = 2.69746994972229;
        c2c = 0.46352824568748474;
        c3a = 2.260690212249756;
        c3b = 1.0386223793029785;
        c3c = 1.1139812469482422;
        L = 28.704402923583984;
        k = 5.6757707595825195;
        x0 = -2.7465617656707764;
        L1 = 0.03422816842794418;
        k1 = 9.087733268737793;
        x_01 = 1.646586537361145;
        L2 = 0.0360025055706501;
        k2 = 13.899486541748047;
        x_02 = -1.636985421180725;
        L3 = 0.07758871465921402;
        k3 = 6.658352375030518;
        x_03 = 1.6404253244400024;

        b1 = 0.0001*(c1a*L4ESeg+c1b*L4ISeg+c1c*L6Ind);
        b2 = 0.0001*(c2a*L4ESeg+c2b*L4ISeg+c2c*L6Ind);
        b3 = 0.0001*(c3a*L4ESeg+c3b*L4ISeg+c3c*L6Ind);

        b1 = L1./(1+exp(-k1*(b1-x_01)));
        b2 = L2./(1+exp(-k2*(b2-x_02)));
        b3 = L3./(1+exp(-k3*(b3-x_03)));

        h = 0.001*a1*(1+b1).*L4ESeg+0.001*a2*(1+b2).*L4ISeg +0.1*a3*(1+b3).*L6Ind;
        LDEOutLIBy_2_20 =  L./(1+exp(-k*0.1*(h-x0)));

       end

      function LDEOutLIBy_10_50 = I_LGNc1_10_50(L4ESeg,L4ISeg,L6Ind)
        % parameter for 10-50Hz
        a1 = 1.2567967176437378;
        a2 = -2.0847339630126953;
        a3 = 2.8438634872436523;
        c1a = 1.9547667503356934;
        c1b = -0.10223595052957535;
        c1c = 0.31979620456695557;
        c2a = 0.6164760589599609;
        c2b = 2.69746994972229;
        c2c = 0.46352824568748474;
        c3a = 2.173144817352295;
        c3b = 0.8922051787376404;
        c3c = 1.2048964500427246;
        L = 78.23750305175781;
        k = 4.4519829750061035;
        x0 = -1.1303282976150513;
        L1 = 0.03653157874941826;
        k1 = 9.096826553344727;
        x_01 = 1.3424837589263916;
        L2 = 0.0378168486058712;
        k2 = 13.899486541748047;
        x_02 = -1.638177514076233;
        L3 = 0.062262631952762604;
        k3 = 6.870923042297363;
        x_03 = 1.3555337190628052;
        b1 = 0.0001*(c1a*L4ESeg+c1b*L4ISeg+c1c*L6Ind);
        b2 = 0.0001*(c2a*L4ESeg+c2b*L4ISeg+c2c*L6Ind);
        b3 = 0.0001*(c3a*L4ESeg+c3b*L4ISeg+c3c*L6Ind);

        b1 = L1./(1+exp(-k1*(b1-x_01)));
        b2 = L2./(1+exp(-k2*(b2-x_02)));
        b3 = L3./(1+exp(-k3*(b3-x_03)));

        h = 0.001*a1*(1+b1).*L4ESeg+0.001*a2*(1+b2).*L4ISeg +0.1*a3*(1+b3).*L6Ind;
        LDEOutLIBy_10_50 =  L./(1+exp(-k*0.1*(h-x0)));

    end


    function LDEOutLIBy_40_100 = I_LGNc1_40_100(L4ESeg,L4ISeg,L6Ind)
        %parameter for 40-100Hz
        a1_1 = 1.2497105598449707;
        a2_1 = -2.1369783878326416;
        a3_1 = 2.8395955562591553;
        c1a_1 = 0.952480673789978;
        c1b_1 = 0.2115955799818039;
        c1c_1 = 0.6714408993721008;
        c2a_1 = 0.6170667409896851;
        c2b_1 = 2.69746994972229;
        c2c_1 = 0.46352824568748474;
        c3a_1 = 1.6640061140060425;
        c3b_1 = 1.4071661233901978;
        c3c_1 = 1.130609154701233;
        L_1 = 126.52735900878906;
        k_1 = 4.514898300170898;
        x0_1 = -0.5162329077720642;
        L1_1 = 0.05974283814430237;
        k1_1 = 9.959728240966797;
        x_01_1 = 0.7407005429267883;
        L2_1 = 0.061291903257369995;
        k2_1 = 13.899486541748047;
        x_02_1 = -1.6503574848175049;
        L3_1 = 0.07718821614980698;
        k3_1 = 8.294986724853516;
        x_03_1 = 1.009885549545288;


        b1_1 = 0.0001*(c1a_1*L4ESeg+c1b_1*L4ISeg+c1c_1*L6Ind);
        b2_1 = 0.0001*(c2a_1*L4ESeg+c2b_1*L4ISeg+c2c_1*L6Ind);
        b3_1 = 0.0001*(c3a_1*L4ESeg+c3b_1*L4ISeg+c3c_1*L6Ind);

        b1_1 = L1_1./(1+exp(-k1_1*(b1_1-x_01_1)));
        b2_1 = L2_1./(1+exp(-k2_1*(b2_1-x_02_1)));
        b3_1 = L3_1./(1+exp(-k3_1*(b3_1-x_03_1)));

        h_1 = 0.001*a1_1*(1+b1_1).*L4ESeg+0.001*a2_1*(1+b2_1).*L4ISeg +0.1*a3_1*(1+b3_1).*L6Ind;
        LDEOutLIBy_40_100 =  L_1./(1+exp(-k_1*0.1*(h_1-x0_1)));
    end

    function LDEOutLIBy_90_160 = I_LGNc1_90_160(L4ESeg,L4ISeg,L6Ind)
           % parameter for 90-160Hz
        a1_1 = 1.2358331680297852;
        a2_1 = -2.1552371978759766;
        a3_1 = 2.820404529571533;
        c1a_1 = 0.3043330907821655;
        c1b_1 = 0.3109455108642578;
        c1c_1 = 10.458990097045898;
        c2a_1 = 0.6199336647987366;
        c2b_1 = 2.6998541355133057;
        c2c_1 = 0.46352824568748474;
        c3a_1 = 0.7396562695503235;
        c3b_1 = 1.3159891366958618;
        c3c_1 = 10.53139591217041;
        L_1 = 179.80624389648438;
        k_1 = 4.5910444259643555;
        x0_1 = 0.683269739151001;
        L1_1 = 0.11333775520324707;
        k1_1 = 14.17318344116211;
        x_01_1 = 0.3612877428531647;
        L2_1 = 0.11951369792222977;
        k2_1 = 13.907645225524902;
        x_02_1 = -1.6939187049865723;
        L3_1 = 0.12188170105218887;
        k3_1 = 14.489945411682129;
        x_03_1 = 0.6155731081962585;


        b1_1 = 0.0001*(c1a_1*L4ESeg+c1b_1*L4ISeg+c1c_1*L6Ind);
        b2_1 = 0.0001*(c2a_1*L4ESeg+c2b_1*L4ISeg+c2c_1*L6Ind);
        b3_1 = 0.0001*(c3a_1*L4ESeg+c3b_1*L4ISeg+c3c_1*L6Ind);

        b1_1 = L1_1./(1+exp(-k1_1*(b1_1-x_01_1)));
        b2_1 = L2_1./(1+exp(-k2_1*(b2_1-x_02_1)));
        b3_1 = L3_1./(1+exp(-k3_1*(b3_1-x_03_1)));

        h_1 = 0.001*a1_1*(1+b1_1).*L4ESeg+0.001*a2_1*(1+b2_1).*L4ISeg +0.1*a3_1*(1+b3_1).*L6Ind;
        LDEOutLIBy_90_160 =  L_1./(1+exp(-k_1*0.1*(h_1-x0_1)));
    end

    function LDEOutLIBy_150_350 = I_LGNc1_150_350(L4ESeg,L4ISeg,L6Ind)
        %parameter for 150-350Hz
        a1_1 = 0.7290038466453552;
        a2_1 = -1.2425990104675293;
        a3_1 = 1.6305961608886719;
        c1a_1 = -0.6129722595214844;
        c1b_1 = 2.002471685409546;
        c1c_1 = 0.4681350290775299;
        c2a_1 = 0.6229845881462097;
        c2b_1 = 2.6998541355133057;
        c2c_1 = 0.46352824568748474;
        c3a_1 = -1.1220356225967407;
        c3b_1 = 2.586432695388794;
        c3c_1 = 0.8974246978759766;
        L_1 = 377.2705383300781;
        k_1 = 2.4571025371551514;
        x0_1 = 3.935105800628662;
        L1_1 = 0.19075772166252136;
        k1_1 = 10.428500175476074;
        x_01_1 = 0.1070801168680191;
        L2_1 = 0.23380036652088165;
        k2_1 = 13.907645225524902;
        x_02_1 = -1.718724012374878;
        L3_1 = 0.20249710977077484;
        k3_1 = 9.155401229858398;
        x_03_1 = 0.3893190026283264;

        b1_1 = 0.0001*(c1a_1*L4ESeg+c1b_1*L4ISeg+c1c_1*L6Ind);
        b2_1 = 0.0001*(c2a_1*L4ESeg+c2b_1*L4ISeg+c2c_1*L6Ind);
        b3_1 = 0.0001*(c3a_1*L4ESeg+c3b_1*L4ISeg+c3c_1*L6Ind);

        b1_1 = L1_1./(1+exp(-k1_1*(b1_1-x_01_1)));
        b2_1 = L2_1./(1+exp(-k2_1*(b2_1-x_02_1)));
        b3_1 = L3_1./(1+exp(-k3_1*(b3_1-x_03_1)));

        h_1 = 0.001*a1_1*(1+b1_1).*L4ESeg+0.001*a2_1*(1+b2_1).*L4ISeg +0.1*a3_1*(1+b3_1).*L6Ind;
        LDEOutLIBy_150_350 =  L_1./(1+exp(-k_1*0.1*(h_1-x0_1)));
 end


