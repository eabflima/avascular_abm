#!/bin/bash
rm *# *~ circle_cell.m &> /dev/null
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#--------------- Generate circle function ---------
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
echo "function circle_cell(t,x,y,r)" >> circle_cell.m
echo "THETA=linspace(0,2*pi,100);" >> circle_cell.m
echo "RHO=ones(1,100)*r;" >> circle_cell.m
echo "[X,Y] = pol2cart(THETA,RHO);" >> circle_cell.m
echo "X=X+x;" >> circle_cell.m
echo "Y=Y+y;" >> circle_cell.m
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#--------------- Tumor cells ----------------------
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#################### Dead tumor cells #########################
echo "if t == 0" >> circle_cell.m
echo "h=fill(X,Y,[0.6 0.6 0.6],'LineWidth',0.01);" >> circle_cell.m
echo "set(h,'facealpha',.75);" >> circle_cell.m
echo "axis square;" >> circle_cell.m
#################### Quiescente tumor cells ###################
echo "elseif t == 1" >> circle_cell.m
echo "h=fill(X,Y,[0.0 0.0 1.0],'LineWidth',0.01);" >> circle_cell.m
echo "set(h,'facealpha',.75);" >> circle_cell.m
echo "axis square;" >> circle_cell.m
#################### Proliferative tumor cells ################
echo "elseif t == 2" >> circle_cell.m
echo "h=fill(X,Y,[0.0 1.0 0.0],'LineWidth',0.01);" >> circle_cell.m
echo "set(h,'facealpha',.75);" >> circle_cell.m
echo "axis square;" >> circle_cell.m
#################### Hypoxic tumor cells ######################
echo "elseif t == 3" >> circle_cell.m
echo "h=fill(X,Y,[0.0 0.0 0.0],'LineWidth',0.01);" >> circle_cell.m
echo "set(h,'facealpha',.75);" >> circle_cell.m
echo "axis square;" >> circle_cell.m
#################### Dead tumor cells #########################
echo "elseif t == 4" >> circle_cell.m
echo "h=fill(X,Y,[1.0 0.0 0.0],'LineWidth',0.01);" >> circle_cell.m
echo "set(h,'facealpha',.75);" >> circle_cell.m
echo "axis square;" >> circle_cell.m
#################### Growing tumor cells ######################
echo "elseif t == 5" >> circle_cell.m
echo "h=fill(X,Y,[1.0 1.0 0.0],'LineWidth',0.01);" >> circle_cell.m
echo "set(h,'facealpha',.75);" >> circle_cell.m
echo "axis square;" >> circle_cell.m
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#--------------- Endothelial cells ----------------
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
echo "elseif t == 7" >> circle_cell.m
echo "h=fill(X,Y,'r','LineWidth',0.01);" >> circle_cell.m
echo "set(h,'facealpha',.75);" >> circle_cell.m
echo "axis square;" >> circle_cell.m
echo "elseif t == 8" >> circle_cell.m
echo "h=fill(X,Y,[0 0.5 0],'LineWidth',0.01);" >> circle_cell.m
echo "set(h,'facealpha',.75);" >> circle_cell.m
echo "axis square;" >> circle_cell.m
echo "elseif t == 9" >> circle_cell.m
echo "h=fill(X,Y,[0 1 1],'LineWidth',0.01);" >> circle_cell.m
echo "set(h,'facealpha',.75);" >> circle_cell.m
echo "elseif t == 15" >> circle_cell.m
echo "h=fill(X,Y,[.1 .1 .1],'LineWidth',0.01);" >> circle_cell.m
echo "set(h,'facealpha',.75);" >> circle_cell.m
echo "axis square;" >> circle_cell.m
echo "elseif t == 12" >> circle_cell.m
echo "h=fill(X,Y,[1 1 0],'LineWidth',0.01);" >> circle_cell.m
echo "set(h,'facealpha',.75);" >> circle_cell.m
echo "axis square;" >> circle_cell.m
echo "elseif t == 13" >> circle_cell.m
echo "h=fill(X,Y,[1 1 0],'LineWidth',0.01);" >> circle_cell.m
echo "set(h,'facealpha',.75);" >> circle_cell.m
echo "axis square;" >> circle_cell.m
echo "elseif t == 14" >> circle_cell.m
echo "h=fill(X,Y,[0 1 0],'LineWidth',0.01);" >> circle_cell.m
echo "set(h,'facealpha',.75);" >> circle_cell.m
echo "axis square;" >> circle_cell.m
echo "elseif t == 10" >> circle_cell.m
echo "h=fill(X,Y,[1 0.843137 0],'LineWidth',0.01);" >> circle_cell.m
echo "set(h,'facealpha',.75);" >> circle_cell.m
echo "axis square;" >> circle_cell.m
echo "else" >> circle_cell.m
echo "ang=0:0.01:2*pi;" >> circle_cell.m
echo "xp=r*cos(ang);" >> circle_cell.m
echo "yp=r*sin(ang);" >> circle_cell.m
echo "plot(x+xp,y+yp,'k');" >> circle_cell.m
echo "end" >> circle_cell.m
echo "end" >> circle_cell.m
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#--------------- Plot files -----------------------
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
rm matlab_figures.m &> /dev/null
echo "clear all;" >> matlab_figures.m
echo "v = VideoWriter('newfile.avi','Uncompressed AVI');" >> matlab_figures.m
echo "v.FrameRate = 5;" >> matlab_figures.m
echo "open(v)" >> matlab_figures.m
for f in saida*.m; do
    file=$(echo ${f} | cut -d. -f1)
    echo "${file};" >> matlab_figures.m
    echo "circle_cell(cells(1,1),cells(1,2),cells(1,3),cells(1,4));" >> matlab_figures.m
    echo "xlim([0 2*cells(1,4)]);" >> matlab_figures.m
    echo "ylim([0 2*cells(1,4)]);" >> matlab_figures.m
    echo "hold on;" >> matlab_figures.m
    echo "s=size(cells);" >> matlab_figures.m
    echo "for c = 2:s(1);" >> matlab_figures.m
    echo "circle_cell(cells(c,1),cells(c,2),cells(c,3),cells(c,4));" >> matlab_figures.m
    echo "end;" >> matlab_figures.m
    echo "print -dpng ${file}.png" >> matlab_figures.m
    echo "A = imread('${file}.png');" >> matlab_figures.m
    echo "writeVideo(v,A)" >> matlab_figures.m
    echo "clf;" >> matlab_figures.m
done
echo "close(v)" >> matlab_figures.m
echo "exit" >> matlab_figures.m
matlab -nodesktop -nosplash -r "matlab_figures"
