function plotDICs(valc,grainsORebsd,grains,grainsORebsd_val,h)
grainsORebsd_val=grainsORebsd_val+2;
if grainsORebsd_val==3%histograms
   histoplot(grainsORebsd,valc,grains,h)
elseif grainsORebsd_val>3%comparisons
     plotComparisons(grainsORebsd,valc,grains,grainsORebsd_val,h)
    
else%plot microstructure
    plotgrainsEBSD(grainsORebsd,valc,grains,h)
end

end

%%
function plotComparisons(grainsORebsd,valc,grains,grainsORebsd_val,h)
grains2 = grains(grainsORebsd.area>50);

switch grainsORebsd_val
    case 4
        %% DSchmid calcs
        ori =grains2('iron fcc').meanOrientation;
        M = zeros(3);M(1,1) = 1;% M(2,2)=-.5;M(3,3)=-.5;%if want to use matrix form
        sigma001 = tensor(M,'name','stress');
        sigmaCS = rotate(sigma001,inv(ori));%inv converts specimen to crystal
        disp('make below more general later')
%         CS=getDBSDpref('CS');CS_=CS{2};
        CS = {crystalSymmetry('m-3m', [3.6599 3.6599 3.6599], 'mineral', 'Iron fcc', 'color', 'light blue'),...
    crystalSymmetry('m-3m', [2.866 2.866 2.866], 'mineral', 'Iron bcc (old)', 'color', [.8 .8 .8])};        
        CS_=CS{1};
        sS = symmetrise(slipSystem.fcc(CS_));
        pl = sS.n; % normal to the slip plane
        b = sS.b; % slip direction in the slip plane
        [tauMax,mActive,nActive,tau,ind] =  calcShearStress(sigmaCS,b,pl); 
        dataEB=abs(tauMax);
        ytext='Schmid Factor';
    case 5%GOS
        dataEB=grains2.GOS;
        ytext='Grain Orientation Spread';
    case 6%grain size
        dataEB=grains2.area;
        ytext='Grain Area';
    case 7
        disp('in progress')
        
    case 8
        disp('soon')
end
     
    
switch valc
    case 1%Exx
        dataH=grains2.Exx;
        xtext='\epsilon_x_x';
        xrange=[0 0.2];

    case 2%Eyy
        dataH=grains2.Eyy;
        xtext='\epsilon_y_y';
        xrange=[-.2  0.1];
        
    case 3%Exy
        dataH=grains2.Exy;
        xtext='\epsilon_x_y'; 
        xrange=[-.2 0.2];
        
    case 4%Eyx
        dataH=grains2.Eyx;
        xtext='\epsilon_y_x';
        xrange=[-.2 0.2];
        
    case 5%EpMax
        dataH=grains2.EpMax;
        xtext='E_p(max)';
        xrange=[0 0.3];
        
    case 6%EpMin
        dataH=grains2.EpMin;
        xtext='E_p(min)';
        xrange=[-.3 0.3];
        
    case 7%Gmax Shear MAx
        dataH=grains2.Gmax;
        xtext='G(max)';
        xrange=[-.3 0.3];
        
    case 8%Theta
        dataH=grains2.Theta;
        xtext='Theta';
        xrange=[-.3 0.3];
        
    case 9%Rotation
        dataH=grains2.Rotation;
        xtext='Rotation';
        xrange=[-.3 0.3];
             
    case 10%Poisson
        dataH=grains2.Poisson;
        xtext='Poissons ratio';
        xrange=[-1 1];
             
end
%h = axes(h);
plot(h,abs(dataH),dataEB,'o');
        xlabel(xtext,'fontangle','italic')
        ylabel(ytext)
%         set(gca,'fontsize',15)
        grid
        hold off
end
%%
function plotgrainsEBSD(grains2,valc,h,pos_,grains)

switch valc
    case 1%Exx
        dataH=grains2.Exx;
        xtext='\epsilon_x_y';
        xrange=[0 0.2];

    case 2%Eyy
        dataH=grains2.Eyy;
        xtext='\epsilon_y_y';
        xrange=[-.2  0.1];
        
    case 3%Exy
        dataH=grains2.Exy;
        xtext='\epsilon_x_y'; 
        xrange=[-.2 0.2];
        
    case 4%Eyx
        dataH=grains2.Eyx;
        xtext='\epsilon_y_x';
        xrange=[-.2 0.2];
        
    case 5%EpMax
        dataH=grains2.EpMax;
        xtext='E_p(max)';
        xrange=[0 0.3];
        
    case 6%EpMin
        dataH=grains2.EpMin;
        xtext='E_p(min)';
        xrange=[-.3 0.3];
        
    case 7%Gmax Shear MAx
        dataH=grains2.Gmax;
        xtext='G(max)';
        xrange=[-.3 0.3];
        
    case 8%Theta
        dataH=grains2.Gmax;
        xtext='Theta';
        xrange=[-.3 0.3];
        
    case 9%Rotation
        dataH=grains2.Rotation;
        xtext='Rotation';
        xrange=[-.3 0.3];
    
    case 10%Poisson
        dataH=grains2.Poisson;
        xtext='Poissons ratio';
        xrange=[-1 1];
end

plot(grains2,dataH,h,'Annotations','Labels')
ylabel(xtext,'fontangle','italic'),xlabel('')
hold on
plot(grains.boundary,'linewidth',2)
hold off
mtexColorMap jet
mtexColorbar
caxis(xrange)
  
end
%%
function histoplot(grainsORebsd,valc,grains,h)
 grains2 = grains(grainsORebsd.area>20);
    switch valc
        case 1
    dataH=grains2.Exx;
    xtext='\epsilon_x_x';
    xrange=-.05:0.01:0.3;
        case 2
    dataH=grains2.Eyy;
    xtext='\epsilon_y_y';
    xrange=-.2:0.01:0.1;
        case 3
    dataH=grains2.Exy;
    xtext='\epsilon_x_y'; 
    xrange=-.2:0.01:0.2;
        case 4
    dataH=grains2.Eyx;
    xtext='\epsilon_y_x';
    xrange=-.2:0.01:0.2;
        case 5%EpMax
    dataH=grains2.EpMax;
    xtext='E_p(max)';
    xrange=0:0.01:0.3;
        case 6%EpMin
    dataH=grains2.EpMin;
    xtext='E_p(min)';
    xrange=-.3:0.01:0.3;
        case 7%Gmax Shear MAx
    dataH=grains2.Gmax;
    xtext='G(max)';
    xrange=-.3:0.01:0.3;
        case 8%Theta
    dataH=grains2.Gmax;
    xtext='Theta';
    xrange=-.3:0.01:0.3;
        case 9%Rotation
    dataH=grains2.Rotation;
    xtext='Rotation';
    xrange=-.3:0.01:0.3;
    case 10%Poisson
        dataH=grains2.Poisson;
        xtext='Poissons ratio';
        xrange=-1:0.1:1;
    end
    %histogram(dataH,xrange,'Normalization','probability');
    hist(h,dataH,xrange)
    xlabel(xtext)
    ylabel('Frequency')
%     set(gca,'fontsize',15)
    grid
    
end