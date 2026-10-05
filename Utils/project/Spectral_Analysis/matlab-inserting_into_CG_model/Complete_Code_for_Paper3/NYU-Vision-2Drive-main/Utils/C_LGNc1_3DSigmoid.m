%% 20250407 Qingyu

function LDEOutLIBy = C_LGNc1_3DSigmoid(L4EUse,L4IUse,L6LibInd)  %#codegen
     

    % this is the rough sigmoid function used to approximate which segment
    % should we use
    a1 = 0.5820584893226624;
    a2 = -0.9147660136222839;
    a3 = 1.9534958600997925;
    c1a = -0.1777229756116867;
    c1b = 0.2791021466255188;
    c1c = -17.00741958618164;
    c2a = 0.5826088190078735;
    c2b = 2.6192309856414795;
    c2c = 0.3780990242958069;
    c3a = -1.10655677318573;
    c3b = 1.416273832321167;
    c3c = -16.743148803710938;
    L = 252.6310272216797;
    k = 4.470812797546387;
    x0 = 6.122119426727295;
    L1 = 0.8673499822616577;
    k1 = 8.222753524780273;
    x_01 = -0.2331233024597168;
    L2 = 0.7344233989715576;
    k2 = 13.78746223449707;
    x_02 = -1.479231595993042;
    L3 = 0.6437055468559265;
    k3 = 4.013402462005615;
    x_03 = -0.11783179640769958;

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
    LDEOutLIBy(idx1) = C_LGNc1_2_20(L4EUse(idx1), L4IUse(idx1), L6LibInd(idx1));
    LDEOutLIBy(idx2) = C_LGNc1_40_100(L4EUse(idx2), L4IUse(idx2), L6LibInd(idx2));
    LDEOutLIBy(idx3) = C_LGNc1_90_160(L4EUse(idx3), L4IUse(idx3), L6LibInd(idx3));
    LDEOutLIBy(idx4) = C_LGNc1_150_250(L4EUse(idx4), L4IUse(idx4), L6LibInd(idx4));
    LDEOutLIBy(idx9) = C_LGNc1_10_50(L4EUse(idx9), L4IUse(idx9), L6LibInd(idx9));


    % Interpolated values

    % here use a weighted interpolation to get the firing rate of the
    % overlapping part.
    LDEOutLIBy(idx5) = ((LDEOutLIBy_rough(idx5)-10)/10) .* C_LGNc1_2_20(L4EUse(idx5), L4IUse(idx5), L6LibInd(idx5)) + ...
                             ((20-LDEOutLIBy_rough(idx5))/10) .* C_LGNc1_10_50(L4EUse(idx5), L4IUse(idx5), L6LibInd(idx5));

    LDEOutLIBy(idx6) = ((LDEOutLIBy_rough(idx6)-40)/10) .* C_LGNc1_10_50(L4EUse(idx6), L4IUse(idx6), L6LibInd(idx6)) + ...
                             ((50-LDEOutLIBy_rough(idx6))/10) .* C_LGNc1_40_100(L4EUse(idx6), L4IUse(idx6), L6LibInd(idx6));

    LDEOutLIBy(idx7) = ((LDEOutLIBy_rough(idx7)-90)/10) .* C_LGNc1_40_100(L4EUse(idx7), L4IUse(idx7), L6LibInd(idx7)) + ...
                             ((100-LDEOutLIBy_rough(idx7))/10) .* C_LGNc1_90_160(L4EUse(idx7), L4IUse(idx7), L6LibInd(idx7));

    LDEOutLIBy(idx8) = ((LDEOutLIBy_rough(idx8)-150)/10) .* C_LGNc1_90_160(L4EUse(idx8), L4IUse(idx8), L6LibInd(idx8)) + ...
                             ((160-LDEOutLIBy_rough(idx8))/10) .* C_LGNc1_150_250(L4EUse(idx8), L4IUse(idx8), L6LibInd(idx8));

end


        % this is for the output firing rate 2-20Hz
       function LDEOutLIBy_2_20 = C_LGNc1_2_20(L4ESeg,L4ISeg,L6Ind)

        a1 = 1.067445993423462;
        a2 = -1.6319366693496704;
        a3 = 3.740851402282715;
        c1a = -0.10239136219024658;
        c1b = 0.13905571401119232;
        c1c = -15.150544166564941;
        c2a = 0.5853922367095947;
        c2b = 2.6235127449035645;
        c2c = 0.3780990242958069;
        c3a = -0.8116466999053955;
        c3b = 0.45561543107032776;
        c3c = -48.63163757324219;
        L = 38.710655212402344;
        k = 5.4058990478515625;
        x0 = 1.977661371231079;
        L1 = 1.6237294673919678;
        k1 = 5.695254325866699;
        x_01 = 0.10110731422901154;
        L2 = 0.4796615242958069;
        k2 = 13.796998977661133;
        x_02 = -1.574683427810669;
        L3 = 0.4492473006248474;
        k3 = 2.757380723953247;
        x_03 = -0.8385544419288635;


        b1 = 0.0001*(c1a*L4ESeg+c1b*L4ISeg+c1c*L6Ind);
        b2 = 0.0001*(c2a*L4ESeg+c2b*L4ISeg+c2c*L6Ind);
        b3 = 0.0001*(c3a*L4ESeg+c3b*L4ISeg+c3c*L6Ind);

        b1 = L1./(1+exp(-k1*(b1-x_01)));
        b2 = L2./(1+exp(-k2*(b2-x_02)));
        b3 = L3./(1+exp(-k3*(b3-x_03)));

        h = 0.001*a1*(1+b1).*L4ESeg+0.001*a2*(1+b2).*L4ISeg +0.1*a3*(1+b3).*L6Ind;
        LDEOutLIBy_2_20 =  L./(1+exp(-k*0.1*(h-x0)));

       end

       % this is for the output firing rate: 10-50Hz
      function LDEOutLIBy_10_50 = C_LGNc1_10_50(L4ESeg,L4ISeg,L6Ind)

        % 5.21 relative mse loss

        a1 = 0.8615175485610962;
        a2 = -1.3681024312973022;
        a3 = 2.9557478427886963;
        c1a = -0.10280900448560715;
        c1b = 0.15140660107135773;
        c1c = -11.41611385345459;
        c2a = 0.5893976092338562;
        c2b = 2.632138729095459;
        c2c = 0.3781377077102661;
        c3a = -0.31536102294921875;
        c3b = 0.23273396492004395;
        c3c = -48.23828125;
        L = 84.6928939819336;
        k = 4.53830623626709;
        x0 = 3.291092872619629;
        L1 = 1.9911925792694092;
        k1 = 3.925595760345459;
        x_01 = 0.17611198127269745;
        L2 = 0.59935063123703;
        k2 = 13.817185401916504;
        x_02 = -1.718879222869873;
        L3 = 0.6644960045814514;
        k3 = 4.846692085266113;
        x_03 = -0.3563142716884613;




        b1 = 0.0001*(c1a*L4ESeg+c1b*L4ISeg+c1c*L6Ind);
        b2 = 0.0001*(c2a*L4ESeg+c2b*L4ISeg+c2c*L6Ind);
        b3 = 0.0001*(c3a*L4ESeg+c3b*L4ISeg+c3c*L6Ind);

        b1 = L1./(1+exp(-k1*(b1-x_01)));
        b2 = L2./(1+exp(-k2*(b2-x_02)));
        b3 = L3./(1+exp(-k3*(b3-x_03)));

        h = 0.001*a1*(1+b1).*L4ESeg+0.001*a2*(1+b2).*L4ISeg +0.1*a3*(1+b3).*L6Ind;
        LDEOutLIBy_10_50 =  L./(1+exp(-k*0.1*(h-x0)));

    end

    % this is for the output firing rate 40-100Hz
    function LDEOutLIBy_40_100 = C_LGNc1_40_100(L4ESeg,L4ISeg,L6Ind)

        % 5.21 using relative mse loss
        a1 = 0.7946094274520874;
        a2 = -1.3366072177886963;
        a3 = 2.744201183319092;
        c1a = -0.12320206314325333;
        c1b = 0.1587771326303482;
        c1c = -34.40204620361328;
        c2a = 0.5873105525970459;
        c2b = 2.623462677001953;
        c2c = 0.37815549969673157;
        c3a = -0.5389071702957153;
        c3b = -0.3120025396347046;
        c3c = -9.718306541442871;
        L = 140.4534454345703;
        k = 4.199314117431641;
        x0 = 3.27001690864563;
        L1 = 0.04447046294808388;
        k1 = 25.763208389282227;
        x_01 = -0.17714834213256836;
        L2 = 0.028292765840888023;
        k2 = 13.796475410461426;
        x_02 = -1.5651638507843018;
        L3 = 0.020764395594596863;
        k3 = 5.208930969238281;
        x_03 = -1.9380253553390503;

        b1_1 = 0.0001*(c1a*L4ESeg+c1b*L4ISeg+c1c*L6Ind);
        b2_1 = 0.0001*(c2a*L4ESeg+c2b*L4ISeg+c2c*L6Ind);
        b3_1 = 0.0001*(c3a*L4ESeg+c3b*L4ISeg+c3c*L6Ind);

        b1_1 = L1./(1+exp(-k1*(b1_1-x_01)));
        b2_1 = L2./(1+exp(-k2*(b2_1-x_02)));
        b3_1 = L3./(1+exp(-k3*(b3_1-x_03)));

        h_1 = 0.001*a1*(1+b1_1).*L4ESeg+0.001*a2*(1+b2_1).*L4ISeg +0.1*a3*(1+b3_1).*L6Ind;
        LDEOutLIBy_40_100 =  L./(1+exp(-k*0.1*(h_1-x0)));
    end

    function LDEOutLIBy_90_160 = C_LGNc1_90_160(L4ESeg,L4ISeg,L6Ind)

        a1_1 = 0.6914945244789124;
        a2_1 = -1.1808913946151733;
        a3_1 = 2.40539813041687;
        c1a_1 = -0.39529287815093994;
        c1b_1 = 1.8488218784332275;
        c1c_1 = -19.218887329101562;
        c2a_1 = 0.5939167141914368;
        c2b_1 = 2.6250550746917725;
        c2c_1 = 0.37849122285842896;
        c3a_1 = -0.6805426478385925;
        c3b_1 = 2.3302719593048096;
        c3c_1 = -16.8259220123291;
        L_1 = 178.55311584472656;
        k_1 = 4.230237007141113;
        x0_1 = 3.7110230922698975;
        L1_1 = 0.1044311672449112;
        k1_1 = 11.309181213378906;
        x_01_1 = 0.2827354669570923;
        L2_1 = 0.12728218734264374;
        k2_1 = 13.802802085876465;
        x_02_1 = -1.5976545810699463;
        L3_1 = 0.10311955213546753;
        k3_1 = 6.310617923736572;
        x_03_1 = 0.7996283173561096;

        b1_1 = 0.0001*(c1a_1*L4ESeg+c1b_1*L4ISeg+c1c_1*L6Ind);
        b2_1 = 0.0001*(c2a_1*L4ESeg+c2b_1*L4ISeg+c2c_1*L6Ind);
        b3_1 = 0.0001*(c3a_1*L4ESeg+c3b_1*L4ISeg+c3c_1*L6Ind);

        b1_1 = L1_1./(1+exp(-k1_1*(b1_1-x_01_1)));
        b2_1 = L2_1./(1+exp(-k2_1*(b2_1-x_02_1)));
        b3_1 = L3_1./(1+exp(-k3_1*(b3_1-x_03_1)));

        h_1 = 0.001*a1_1*(1+b1_1).*L4ESeg+0.001*a2_1*(1+b2_1).*L4ISeg +0.1*a3_1*(1+b3_1).*L6Ind;
        LDEOutLIBy_90_160 =  L_1./(1+exp(-k_1*0.1*(h_1-x0_1)));
    end

    
    function LDEOutLIBy_150_250 = C_LGNc1_150_250(L4ESeg,L4ISeg,L6Ind)

        a1_1 = 0.4685199558734894;
        a2_1 = -0.7631587386131287;
        a3_1 = 1.61995530128479;
        c1a_1 = -0.3154031038284302;
        c1b_1 = 1.287166714668274;
        c1c_1 = -20.431594848632812;
        c2a_1 = 0.5972834825515747;
        c2b_1 = 2.622128963470459;
        c2c_1 = 0.3784196972846985;
        c3a_1 = -0.5399700403213501;
        c3b_1 = 1.3780064582824707;
        c3c_1 = -20.353286743164062;
        L_1 = 269.791015625;
        k_1 = 2.9197447299957275;
        x0_1 = 4.688855171203613;
        L1_1 = 0.1836596131324768;
        k1_1 = 13.03409481048584;
        x_01_1 = 0.06765943765640259;
        L2_1 = 0.2637574374675751;
        k2_1 = 13.797860145568848;
        x_02_1 = -1.5457819700241089;
        L3_1 = 0.171334907412529;
        k3_1 = 9.458685874938965;
        x_03_1 = 0.21971768140792847;

        b1_1 = 0.0001*(c1a_1*L4ESeg+c1b_1*L4ISeg+c1c_1*L6Ind);
        b2_1 = 0.0001*(c2a_1*L4ESeg+c2b_1*L4ISeg+c2c_1*L6Ind);
        b3_1 = 0.0001*(c3a_1*L4ESeg+c3b_1*L4ISeg+c3c_1*L6Ind);

        b1_1 = L1_1./(1+exp(-k1_1*(b1_1-x_01_1)));
        b2_1 = L2_1./(1+exp(-k2_1*(b2_1-x_02_1)));
        b3_1 = L3_1./(1+exp(-k3_1*(b3_1-x_03_1)));

        h_1 = 0.001*a1_1*(1+b1_1).*L4ESeg+0.001*a2_1*(1+b2_1).*L4ISeg +0.1*a3_1*(1+b3_1).*L6Ind;
        LDEOutLIBy_150_250 =  L_1./(1+exp(-k_1*0.1*(h_1-x0_1)));
 end

   
