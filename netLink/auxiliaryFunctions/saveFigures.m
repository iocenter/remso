hfigs = get(0, 'children');                            %Get list of figures

for m = 1:length(hfigs)
    figure(hfigs(m))                                   %Bring Figure to foreground
    
    saveas(hfigs(m), [strcat('figs/', num2str(hfigs(m).Number)), '.fig']) %Matlab .FIG file
    saveas(hfigs(m), [strcat('figs/', num2str(hfigs(m).Number)), '.svg']) % SVG format
    eval(['print -depsc2 -r300 ' strcat('figs/', num2str(hfigs(m).Number))])   %Enhanced Postscript (Level 2 color) (Best for LaTeX documents)    
end