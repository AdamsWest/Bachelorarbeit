figure
hold on
for K_V_max = 200:200:2000
    m_Mot4 = 0;
    y14 = 0;
    for i = 1:length(m_Mot)
        if K_V(i) < K_V_max && K_V(i) > K_V_max - 200
            m_Mot4(i) = m_Mot(i);
            y14(i) = y1(i);
        else
            m_Mot4(i) = 0;
            y14(i) = 0;
        end
    end
    plot(m_Mot4,y14,'o','Color',[1 1 1]*K_V_max/2000)
    hold on
end