for i = 1:6


    p = pie(runtimes198(10,[3,4,6])*1000);
    pText = findobj(p,'Type','text');
    percentValues = get(pText,'String');
    txt = {'I/O Time: ';'Iteration Time: ';'Gather time: '};
    combinedtxt = strcat(txt,percentValues);
    pText(1).String = combinedtxt(1);
    pText(2).String = combinedtxt(2);
    pText(3).String = combinedtxt(3);

    p(2).FontSize = 14;
    p(4).FontSize = 14;
    p(6).FontSize = 14;
    
    cols = [116 201 229; 179 219 210; 60 186 198];
    cols = cols ./ 255;
    p(1).FaceColor = cols(1,:);
    p(3).FaceColor = cols(2,:);
    p(3).FaceColor = cols(3,:);
    
    % Placeholder debug
    
end