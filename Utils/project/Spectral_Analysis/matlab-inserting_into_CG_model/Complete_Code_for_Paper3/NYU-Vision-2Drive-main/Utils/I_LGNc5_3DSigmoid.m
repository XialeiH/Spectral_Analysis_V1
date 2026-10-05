%% 20250407 Qingyu (data correct)

function LDEOutLIBy = I_LGNc5_3DSigmoid(L4EUse,L4IUse,L6LibInd)  %#codegen
  
    % first we predict which segment(s) we should use by a rough estimation

    a1 = 0.8341930508613586;
    a2 = -1.4317063093185425;
    a3 = 1.8338639736175537;
    L = 358.79840087890625;
    k = 3.181288480758667;
    x0 = 5.277807235717773;
            
    c1a = -0.24781513214111328;
    c1b = 0.5277940630912781;
    c1c = -26.88960838317871;
    c2a = 0.6228991746902466;
    c2b = 2.701118230819702;
    c2c = 0.4637240469455719;
    c3a = -0.3395082354545593;
    c3b = 0.5909067392349243;
    c3c = -27.56024742126465;
    L1 = 0.22450560331344604;
    k1 = 23.348894119262695;
    x_01 = -0.12918516993522644;
    L2 = 0.24363087117671967;
    k2 = 13.910262107849121;
    x_02 = -1.7132883071899414;
    L3 = 0.2637702226638794;
    k3 = 23.985416412353516;
    x_03 = -0.06817109882831573;

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
    LDEOutLIBy(idx1) = I_LGNc5_2_20(L4EUse(idx1), L4IUse(idx1), L6LibInd(idx1));
    LDEOutLIBy(idx2) = I_LGNc5_40_100(L4EUse(idx2), L4IUse(idx2), L6LibInd(idx2));
    LDEOutLIBy(idx3) = I_LGNc5_90_160(L4EUse(idx3), L4IUse(idx3), L6LibInd(idx3));
    LDEOutLIBy(idx4) = I_LGNc5_150_350(L4EUse(idx4), L4IUse(idx4), L6LibInd(idx4));
    LDEOutLIBy(idx9) = I_LGNc5_10_50(L4EUse(idx9), L4IUse(idx9), L6LibInd(idx9));

     % Interpolated values
    LDEOutLIBy(idx5) = ((LDEOutLIBy_rough(idx5)-10)/10) .* I_LGNc5_2_20(L4EUse(idx5), L4IUse(idx5), L6LibInd(idx5)) + ...
                             ((20-LDEOutLIBy_rough(idx5))/10) .* I_LGNc5_10_50(L4EUse(idx5), L4IUse(idx5), L6LibInd(idx5));

    LDEOutLIBy(idx6) = ((LDEOutLIBy_rough(idx6)-40)/10) .* I_LGNc5_10_50(L4EUse(idx6), L4IUse(idx6), L6LibInd(idx6)) + ...
                             ((50-LDEOutLIBy_rough(idx6))/10) .* I_LGNc5_40_100(L4EUse(idx6), L4IUse(idx6), L6LibInd(idx6));

    LDEOutLIBy(idx7) = ((LDEOutLIBy_rough(idx7)-90)/10) .* I_LGNc5_40_100(L4EUse(idx7), L4IUse(idx7), L6LibInd(idx7)) + ...
                             ((100-LDEOutLIBy_rough(idx7))/10) .* I_LGNc5_90_160(L4EUse(idx7), L4IUse(idx7), L6LibInd(idx7));

    LDEOutLIBy(idx8) = ((LDEOutLIBy_rough(idx8)-150)/10) .* I_LGNc5_90_160(L4EUse(idx8), L4IUse(idx8), L6LibInd(idx8)) + ...
                             ((160-LDEOutLIBy_rough(idx8))/10) .* I_LGNc5_150_350(L4EUse(idx8), L4IUse(idx8), L6LibInd(idx8));

end

    function LDEOutLIBy_2_20 = I_LGNc5_2_20(L4ESeg,L4ISeg,L6Ind)
        % parameter for 0-20Hz
        a1 = 1.384969711303711;
        a2 = -2.204441547393799;
        a3 = 3.160212993621826;
        c1a = 2.347598075866699;
        c1b = 0.19028203189373016;
        c1c = 0.7332016825675964;
        c2a = 0.6164997816085815;
        c2b = 2.69746994972229;
        c2c = 0.46352824568748474;
        c3a = -0.2402435541152954;
        c3b = 0.2693817913532257;
        c3c = -36.163883209228516;
        L = 28.624727249145508;
        k = 5.2093892097473145;
        x0 = -1.214694857597351;
        L1 = 0.21182674169540405;
        k1 = 5.292797565460205;
        x_01 = 1.0927011966705322;
        L2 = 0.19623348116874695;
        k2 = 13.899486541748047;
        x_02 = -1.644126057624817;
        L3 = 0.2588055729866028;
        k3 = 29.938589096069336;
        x_03 = -0.2255963683128357;

        b1 = 0.0001*(c1a*L4ESeg+c1b*L4ISeg+c1c*L6Ind);
        b2 = 0.0001*(c2a*L4ESeg+c2b*L4ISeg+c2c*L6Ind);
        b3 = 0.0001*(c3a*L4ESeg+c3b*L4ISeg+c3c*L6Ind);

        b1 = L1./(1+exp(-k1*(b1-x_01)));
        b2 = L2./(1+exp(-k2*(b2-x_02)));
        b3 = L3./(1+exp(-k3*(b3-x_03)));

        h = 0.001*a1*(1+b1).*L4ESeg+0.001*a2*(1+b2).*L4ISeg +0.1*a3*(1+b3).*L6Ind;
        LDEOutLIBy_2_20 =  L./(1+exp(-k*0.1*(h-x0)));

       end

      function LDEOutLIBy_10_50 = I_LGNc5_10_50(L4ESeg,L4ISeg,L6Ind)
        % parameter for 10-50Hz
        a1 = 1.3338099718093872;
        a2 = -2.205650568008423;
        a3 = 3.020165205001831;
        c1a = 1.7926210165023804;
        c1b = 0.5422835946083069;
        c1c = 0.6585255861282349;
        c2a = 0.6170958280563354;
        c2b = 2.69746994972229;
        c2c = 0.46352824568748474;
        c3a = -0.2875407636165619;
        c3b = 0.37620577216148376;
        c3c = -36.673431396484375;
        L = 78.09541320800781;
        k = 4.218593597412109;
        x0 = 0.570160984992981;
        L1 = 0.17700161039829254;
        k1 = 4.06898307800293;
        x_01 = 1.2156531810760498;
        L2 = 0.1688620001077652;
        k2 = 13.899486541748047;
        x_02 = -1.6536628007888794;
        L3 = 0.19656233489513397;
        k3 = 30.26332664489746;
        x_03 = -0.19047988951206207;
        b1 = 0.0001*(c1a*L4ESeg+c1b*L4ISeg+c1c*L6Ind);
        b2 = 0.0001*(c2a*L4ESeg+c2b*L4ISeg+c2c*L6Ind);
        b3 = 0.0001*(c3a*L4ESeg+c3b*L4ISeg+c3c*L6Ind);

        b1 = L1./(1+exp(-k1*(b1-x_01)));
        b2 = L2./(1+exp(-k2*(b2-x_02)));
        b3 = L3./(1+exp(-k3*(b3-x_03)));

        h = 0.001*a1*(1+b1).*L4ESeg+0.001*a2*(1+b2).*L4ISeg +0.1*a3*(1+b3).*L6Ind;
        LDEOutLIBy_10_50 =  L./(1+exp(-k*0.1*(h-x0)));

    end
    function LDEOutLIBy_40_100 = I_LGNc5_40_100(L4ESeg,L4ISeg,L6Ind)
        %parameter for 40-100Hz

        a1_1 = 1.304731011390686;
        a2_1 = -2.2269785404205322;
        a3_1 = 2.9620110988616943;
        c1a_1 = 0.9269095659255981;
        c1b_1 = 0.8119890689849854;
        c1c_1 = 0.27975067496299744;
        c2a_1 = 0.6179384589195251;
        c2b_1 = 2.6998541355133057;
        c2c_1 = 0.46352824568748474;
        c3a_1 = -0.34343627095222473;
        c3b_1 = 0.4541914463043213;
        c3c_1 = -35.89155197143555;
        L_1 = 133.87088012695312;
        k_1 = 4.110398769378662;
        x0_1 = 1.4260530471801758;
        L1_1 = 0.10678518563508987;
        k1_1 = 5.030115127563477;
        x_01_1 = 1.027201771736145;
        L2_1 = 0.10551595687866211;
        k2_1 = 13.899953842163086;
        x_02_1 = -1.6710052490234375;
        L3_1 = 0.10645925998687744;
        k3_1 = 29.05759620666504;
        x_03_1 = -0.2786136269569397;


        b1_1 = 0.0001*(c1a_1*L4ESeg+c1b_1*L4ISeg+c1c_1*L6Ind);
        b2_1 = 0.0001*(c2a_1*L4ESeg+c2b_1*L4ISeg+c2c_1*L6Ind);
        b3_1 = 0.0001*(c3a_1*L4ESeg+c3b_1*L4ISeg+c3c_1*L6Ind);

        b1_1 = L1_1./(1+exp(-k1_1*(b1_1-x_01_1)));
        b2_1 = L2_1./(1+exp(-k2_1*(b2_1-x_02_1)));
        b3_1 = L3_1./(1+exp(-k3_1*(b3_1-x_03_1)));

        h_1 = 0.001*a1_1*(1+b1_1).*L4ESeg+0.001*a2_1*(1+b2_1).*L4ISeg +0.1*a3_1*(1+b3_1).*L6Ind;
        LDEOutLIBy_40_100 =  L_1./(1+exp(-k_1*0.1*(h_1-x0_1)));
    end

    function LDEOutLIBy_90_160 = I_LGNc5_90_160(L4ESeg,L4ISeg,L6Ind)
           % parameter for 90-160Hz
        a1_1 = 1.2818405628204346;
        a2_1 = -2.2324347496032715;
        a3_1 = 2.921194553375244;
        c1a_1 = 0.3395068347454071;
        c1b_1 = 0.8570632934570312;
        c1c_1 = 9.444347381591797;
        c2a_1 = 0.6199607253074646;
        c2b_1 = 2.7014100551605225;
        c2c_1 = 0.46352824568748474;
        c3a_1 = -0.36326783895492554;
        c3b_1 = 0.47401347756385803;
        c3c_1 = -42.53866958618164;
        L_1 = 179.77674865722656;
        k_1 = 4.510212421417236;
        x0_1 = 2.334826707839966;
        L1_1 = 0.06778999418020248;
        k1_1 = 9.603436470031738;
        x_01_1 = 0.7094544172286987;
        L2_1 = 0.06973715126514435;
        k2_1 = 13.900339126586914;
        x_02_1 = -1.6965258121490479;
        L3_1 = 0.059314001351594925;
        k3_1 = 37.33075714111328;
        x_03_1 = -0.39231711626052856;


        b1_1 = 0.0001*(c1a_1*L4ESeg+c1b_1*L4ISeg+c1c_1*L6Ind);
        b2_1 = 0.0001*(c2a_1*L4ESeg+c2b_1*L4ISeg+c2c_1*L6Ind);
        b3_1 = 0.0001*(c3a_1*L4ESeg+c3b_1*L4ISeg+c3c_1*L6Ind);

        b1_1 = L1_1./(1+exp(-k1_1*(b1_1-x_01_1)));
        b2_1 = L2_1./(1+exp(-k2_1*(b2_1-x_02_1)));
        b3_1 = L3_1./(1+exp(-k3_1*(b3_1-x_03_1)));

        h_1 = 0.001*a1_1*(1+b1_1).*L4ESeg+0.001*a2_1*(1+b2_1).*L4ISeg +0.1*a3_1*(1+b3_1).*L6Ind;
        LDEOutLIBy_90_160 =  L_1./(1+exp(-k_1*0.1*(h_1-x0_1)));
    end

    function LDEOutLIBy_150_350 = I_LGNc5_150_350(L4ESeg,L4ISeg,L6Ind)
        % parameter for 150-350Hz
        a1_1 = 0.7742679119110107;
        a2_1 = -1.323077917098999;
        a3_1 = 1.7382783889770508;
        c1a_1 = -0.7668265700340271;
        c1b_1 = 3.0358164310455322;
        c1c_1 = -0.5870451331138611;
        c2a_1 = 0.6235365271568298;
        c2b_1 = 2.7014100551605225;
        c2c_1 = 0.46352824568748474;
        c3a_1 = -0.2178926169872284;
        c3b_1 = 0.4876348078250885;
        c3c_1 = -44.645442962646484;
        L_1 = 377.4723205566406;
        k_1 = 2.3570971488952637;
        x0_1 = 5.080848693847656;
        L1_1 = 0.18927721679210663;
        k1_1 = 5.098660945892334;
        x_01_1 = 0.3351143002510071;
        L2_1 = 0.24037592113018036;
        k2_1 = 13.900339126586914;
        x_02_1 = -1.721658706665039;
        L3_1 = 0.21629641950130463;
        k3_1 = 40.4277458190918;
        x_03_1 = -0.08736801147460938;

        b1_1 = 0.0001*(c1a_1*L4ESeg+c1b_1*L4ISeg+c1c_1*L6Ind);
        b2_1 = 0.0001*(c2a_1*L4ESeg+c2b_1*L4ISeg+c2c_1*L6Ind);
        b3_1 = 0.0001*(c3a_1*L4ESeg+c3b_1*L4ISeg+c3c_1*L6Ind);

        b1_1 = L1_1./(1+exp(-k1_1*(b1_1-x_01_1)));
        b2_1 = L2_1./(1+exp(-k2_1*(b2_1-x_02_1)));
        b3_1 = L3_1./(1+exp(-k3_1*(b3_1-x_03_1)));

        h_1 = 0.001*a1_1*(1+b1_1).*L4ESeg+0.001*a2_1*(1+b2_1).*L4ISeg +0.1*a3_1*(1+b3_1).*L6Ind;
        LDEOutLIBy_150_350 =  L_1./(1+exp(-k_1*0.1*(h_1-x0_1)));
 end


