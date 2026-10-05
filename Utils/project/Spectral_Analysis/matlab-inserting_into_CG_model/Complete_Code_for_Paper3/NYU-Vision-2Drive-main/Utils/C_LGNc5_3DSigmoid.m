%% 20250407 Qingyu (correct data, 20250409update)

function LDEOutLIBy = C_LGNc5_3DSigmoid(L4EUse,L4IUse,L6LibInd)  %#codegen
     

    % this is the rough sigmoid function used to approximate which segment
    % we should use
    a1 = 0.6106534004211426;
    a2 = -0.8911330103874207;
    a3 = 1.6762933731079102;
    c1a = -0.23054710030555725;
    c1b = 0.3523716926574707;
    c1c = -13.081350326538086;
    c2a = 0.5897137522697449;
    c2b = 2.6364216804504395;
    c2c = 0.37881481647491455;
    c3a = -0.5599619150161743;
    c3b = 0.7008721828460693;
    c3c = -20.329145431518555;
    L = 234.2545166015625;
    k = 4.909726142883301;
    x0 = 5.977545738220215;
    L1 = 0.7097715139389038;
    k1 = 6.44271183013916;
    x_01 = -0.28981977701187134;
    L2 = 0.6944610476493835;
    k2 = 13.8255033493042;
    x_02 = -1.715758204460144;
    L3 = 0.8553881645202637;
    k3 = 6.644332408905029;
    x_03 = -0.1342685967683792;

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
    LDEOutLIBy(idx1) = C_LGNc5_2_20(L4EUse(idx1), L4IUse(idx1), L6LibInd(idx1));
    LDEOutLIBy(idx2) = C_LGNc5_40_100(L4EUse(idx2), L4IUse(idx2), L6LibInd(idx2));
    LDEOutLIBy(idx3) = C_LGNc5_90_160(L4EUse(idx3), L4IUse(idx3), L6LibInd(idx3));
    LDEOutLIBy(idx4) = C_LGNc5_150_250(L4EUse(idx4), L4IUse(idx4), L6LibInd(idx4));
    LDEOutLIBy(idx9) = C_LGNc5_10_50(L4EUse(idx9), L4IUse(idx9), L6LibInd(idx9));


 % Interpolated values

% here use a weighted interpolation to get the firing rate of the
% overlapping part.
    LDEOutLIBy(idx5) = ((LDEOutLIBy_rough(idx5)-10)/10) .* C_LGNc5_2_20(L4EUse(idx5), L4IUse(idx5), L6LibInd(idx5)) + ...
                             ((20-LDEOutLIBy_rough(idx5))/10) .* C_LGNc5_10_50(L4EUse(idx5), L4IUse(idx5), L6LibInd(idx5));

    LDEOutLIBy(idx6) = ((LDEOutLIBy_rough(idx6)-40)/10) .* C_LGNc5_10_50(L4EUse(idx6), L4IUse(idx6), L6LibInd(idx6)) + ...
                             ((50-LDEOutLIBy_rough(idx6))/10) .* C_LGNc5_40_100(L4EUse(idx6), L4IUse(idx6), L6LibInd(idx6));

    LDEOutLIBy(idx7) = ((LDEOutLIBy_rough(idx7)-90)/10) .* C_LGNc5_40_100(L4EUse(idx7), L4IUse(idx7), L6LibInd(idx7)) + ...
                             ((100-LDEOutLIBy_rough(idx7))/10) .* C_LGNc5_90_160(L4EUse(idx7), L4IUse(idx7), L6LibInd(idx7));

    LDEOutLIBy(idx8) = ((LDEOutLIBy_rough(idx8)-150)/10) .* C_LGNc5_90_160(L4EUse(idx8), L4IUse(idx8), L6LibInd(idx8)) + ...
                             ((160-LDEOutLIBy_rough(idx8))/10) .* C_LGNc5_150_250(L4EUse(idx8), L4IUse(idx8), L6LibInd(idx8));

end



     function LDEOutLIBy_2_20 = C_LGNc5_2_20(L4ESeg,L4ISeg,L6Ind)
      
        a1 = 1.0180405378341675;
        a2 = -1.555922508239746;
        a3 = 3.64280104637146;
        c1a = -0.09297386556863785;
        c1b = 0.1252441257238388;
        c1c = -15.658687591552734;
        c2a = 0.5893293023109436;
        c2b = 2.632004976272583;
        c2c = 0.3780990242958069;
        c3a = -1.0518933534622192;
        c3b = 0.32011905312538147;
        c3c = -34.6608772277832;
        L = 33.90888977050781;
        k = 5.893667697906494;
        x0 = 2.407787561416626;
        L1 = 1.4376496076583862;
        k1 = 6.544376373291016;
        x_01 = 0.0821276307106018;
        L2 = 0.4132365584373474;
        k2 = 13.814549446105957;
        x_02 = -1.6902718544006348;
        L3 = 0.358633816242218;
        k3 = 1.8057782649993896;
        x_03 = -1.5424884557724;
        

        b1 = 0.0001*(c1a*L4ESeg+c1b*L4ISeg+c1c*L6Ind);
        b2 = 0.0001*(c2a*L4ESeg+c2b*L4ISeg+c2c*L6Ind);
        b3 = 0.0001*(c3a*L4ESeg+c3b*L4ISeg+c3c*L6Ind);

        b1 = L1./(1+exp(-k1*(b1-x_01)));
        b2 = L2./(1+exp(-k2*(b2-x_02)));
        b3 = L3./(1+exp(-k3*(b3-x_03)));

        h = 0.001*a1*(1+b1).*L4ESeg+0.001*a2*(1+b2).*L4ISeg +0.1*a3*(1+b3).*L6Ind;
        LDEOutLIBy_2_20 =  L./(1+exp(-k*0.1*(h-x0)));

       end

      function LDEOutLIBy_10_50 = C_LGNc5_10_50(L4ESeg,L4ISeg,L6Ind)
      

        a1 = 0.8831879496574402;
        a2 = -1.4129109382629395;
        a3 = 3.1161892414093018;
        c1a = -0.08651598542928696;
        c1b = 0.12443941086530685;
        c1c = -18.66737937927246;
        c2a = 0.5921487212181091;
        c2b = 2.6311733722686768;
        c2c = 0.3780990242958069;
        c3a = -1.025649905204773;
        c3b = -0.15954944491386414;
        c3c = -26.538562774658203;
        L = 78.91669464111328;
        k = 5.645089626312256;
        x0 = 3.328794479370117;
        L1 = 1.1575614213943481;
        k1 = 7.873226642608643;
        x_01 = 0.11685433238744736;
        L2 = 0.23791196942329407;
        k2 = 13.814475059509277;
        x_02 = -1.681871771812439;
        L3 = 0.2414105236530304;
        k3 = 1.0845931768417358;
        x_03 = -2.0774779319763184;


        b1 = 0.0001*(c1a*L4ESeg+c1b*L4ISeg+c1c*L6Ind);
        b2 = 0.0001*(c2a*L4ESeg+c2b*L4ISeg+c2c*L6Ind);
        b3 = 0.0001*(c3a*L4ESeg+c3b*L4ISeg+c3c*L6Ind);

        b1 = L1./(1+exp(-k1*(b1-x_01)));
        b2 = L2./(1+exp(-k2*(b2-x_02)));
        b3 = L3./(1+exp(-k3*(b3-x_03)));

        h = 0.001*a1*(1+b1).*L4ESeg+0.001*a2*(1+b2).*L4ISeg +0.1*a3*(1+b3).*L6Ind;
        LDEOutLIBy_10_50 =  L./(1+exp(-k*0.1*(h-x0)));

    end


    function LDEOutLIBy_40_100 = C_LGNc5_40_100(L4ESeg,L4ISeg,L6Ind)

        a1_1 = 0.95555579662323;
        a2_1 = -1.4942771196365356;
        a3_1 = 2.9723029136657715;
        c1a_1 = -0.060442421585321426;
        c1b_1 = 0.09377294033765793;
        c1c_1 = -23.825483322143555;
        c2a_1 = 0.612342357635498;
        c2b_1 = 2.7010326385498047;
        c2c_1 = 0.37836408615112305;
        c3a_1 = 5.626166820526123;
        c3b_1 = 0.4467492401599884;
        c3c_1 = -19.8721981048584;
        L_1 = 123.07035827636719;
        k_1 = 4.379006385803223;
        x0_1 = 5.084381103515625;
        L1_1 = 0.6466155052185059;
        k1_1 = 12.6839599609375;
        x_01_1 = 0.0010086840484291315;
        L2_1 = 0.31195151805877686;
        k2_1 = 13.9404878616333;
        x_02_1 = -2.0276966094970703;
        L3_1 = 0.41925880312919617;
        k3_1 = 13.110725402832031;
        x_03_1 = 0.13788868486881256;

        b1_1 = 0.0001*(c1a_1*L4ESeg+c1b_1*L4ISeg+c1c_1*L6Ind);
        b2_1 = 0.0001*(c2a_1*L4ESeg+c2b_1*L4ISeg+c2c_1*L6Ind);
        b3_1 = 0.0001*(c3a_1*L4ESeg+c3b_1*L4ISeg+c3c_1*L6Ind);

        b1_1 = L1_1./(1+exp(-k1_1*(b1_1-x_01_1)));
        b2_1 = L2_1./(1+exp(-k2_1*(b2_1-x_02_1)));
        b3_1 = L3_1./(1+exp(-k3_1*(b3_1-x_03_1)));

        h_1 = 0.001*a1_1*(1+b1_1).*L4ESeg+0.001*a2_1*(1+b2_1).*L4ISeg +0.1*a3_1*(1+b3_1).*L6Ind;
        LDEOutLIBy_40_100 =  L_1./(1+exp(-k_1*0.1*(h_1-x0_1)));
    end

    function LDEOutLIBy_90_160 = C_LGNc5_90_160(L4ESeg,L4ISeg,L6Ind)
    
        a1_1 = 0.6682639718055725;
        a2_1 = -1.1407105922698975;
        a3_1 = 2.3222317695617676;
        c1a_1 = -0.383513480424881;
        c1b_1 = 1.8392771482467651;
        c1c_1 = -19.03392219543457;
        c2a_1 = 0.5941360592842102;
        c2b_1 = 2.6249122619628906;
        c2c_1 = 0.37849754095077515;
        c3a_1 = -0.7063209414482117;
        c3b_1 = 2.3735389709472656;
        c3c_1 = -16.371870040893555;
        L_1 = 178.6609649658203;
        k_1 = 4.377488136291504;
        x0_1 = 3.965269088745117;
        L1_1 = 0.10105859488248825;
        k1_1 = 11.35363483428955;
        x_01_1 = 0.28017598390579224;
        L2_1 = 0.12497309595346451;
        k2_1 = 13.802659034729004;
        x_02_1 = -1.59611976146698;
        L3_1 = 0.09917835891246796;
        k3_1 = 6.180002689361572;
        x_03_1 = 0.7886167764663696;

        b1_1 = 0.0001*(c1a_1*L4ESeg+c1b_1*L4ISeg+c1c_1*L6Ind);
        b2_1 = 0.0001*(c2a_1*L4ESeg+c2b_1*L4ISeg+c2c_1*L6Ind);
        b3_1 = 0.0001*(c3a_1*L4ESeg+c3b_1*L4ISeg+c3c_1*L6Ind);

        b1_1 = L1_1./(1+exp(-k1_1*(b1_1-x_01_1)));
        b2_1 = L2_1./(1+exp(-k2_1*(b2_1-x_02_1)));
        b3_1 = L3_1./(1+exp(-k3_1*(b3_1-x_03_1)));

        h_1 = 0.001*a1_1*(1+b1_1).*L4ESeg+0.001*a2_1*(1+b2_1).*L4ISeg +0.1*a3_1*(1+b3_1).*L6Ind;
        LDEOutLIBy_90_160 =  L_1./(1+exp(-k_1*0.1*(h_1-x0_1)));
    end

    function LDEOutLIBy_150_250 = C_LGNc5_150_250(L4ESeg,L4ISeg,L6Ind)

        a1_1 = 0.47295302152633667;
        a2_1 = -0.7710962891578674;
        a3_1 = 1.6340621709823608;
        c1a_1 = -0.32391610741615295;
        c1b_1 = 1.3314367532730103;
        c1c_1 = -20.15659523010254;
        c2a_1 = 0.5970264673233032;
        c2b_1 = 2.622001886367798;
        c2c_1 = 0.37837958335876465;
        c3a_1 = -0.547349214553833;
        c3b_1 = 1.3995949029922485;
        c3c_1 = -20.356365203857422;
        L_1 = 268.79302978515625;
        k_1 = 2.912888765335083;
        x0_1 = 4.95598840713501;
        L1_1 = 0.18132072687149048;
        k1_1 = 12.774676322937012;
        x_01_1 = 0.06931983679533005;
        L2_1 = 0.26252853870391846;
        k2_1 = 13.797679901123047;
        x_02_1 = -1.5408272743225098;
        L3_1 = 0.16747964918613434;
        k3_1 = 9.310758590698242;
        x_03_1 = 0.2144562155008316;

        b1_1 = 0.0001*(c1a_1*L4ESeg+c1b_1*L4ISeg+c1c_1*L6Ind);
        b2_1 = 0.0001*(c2a_1*L4ESeg+c2b_1*L4ISeg+c2c_1*L6Ind);
        b3_1 = 0.0001*(c3a_1*L4ESeg+c3b_1*L4ISeg+c3c_1*L6Ind);

        b1_1 = L1_1./(1+exp(-k1_1*(b1_1-x_01_1)));
        b2_1 = L2_1./(1+exp(-k2_1*(b2_1-x_02_1)));
        b3_1 = L3_1./(1+exp(-k3_1*(b3_1-x_03_1)));

        h_1 = 0.001*a1_1*(1+b1_1).*L4ESeg+0.001*a2_1*(1+b2_1).*L4ISeg +0.1*a3_1*(1+b3_1).*L6Ind;
        LDEOutLIBy_150_250 =  L_1./(1+exp(-k_1*0.1*(h_1-x0_1)));
 end


