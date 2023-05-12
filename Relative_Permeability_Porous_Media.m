clc
clear
%%This code is developed for two-phase flow in 2D porous media based on
%%Lattice Boltzman Phase-Field Free Energy model.>>>>>For Demonstartion
%%of Relative Permeability in Porous Media....By Ali Sadeghkhani
%% A: Initialization  
% Describing the initial parameters
w=[4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
cx=[0 1 0 -1 0 1 -1 -1 1];
cy=[0 0 1 0 -1 1 1 -1 -1];
opp=[1 4 5 2 3 8 9 6 7];
maxT=200000;
tplot=100;
conv=1*10^-9;

%saturation and wettability of the wet phase
sw=0.8;
contact_angle=30;

snw=1-sw;
% Prameters data
xmin=1; xmax=200; ymin=1; ymax=200;
AA=zeros(xmax,ymax);
BB=zeros(xmax,ymax);

% Relaxation Time (Viscosity Definition) related to wet fluid
tau1=0.7;  
rho1=1;
dynvisc1=((tau1-0.5)/3);
matdynvisc1=dynvisc1*ones(xmax,ymax);

% Relaxation Time (Viscosity Definition) related to non-wet-fluid
tau2=1.7;  
rho2=1;
dynvisc2=((tau2-0.5)/3);
matdynvisc2=dynvisc2*ones(xmax,ymax);

% Fluid flow obstacle definition
obst=imread('C:\Users\SONY\Desktop\LBM Study\Porous Media\200_200_Periodic.bmp');
obst=im2bw(obst);
imshow(obst);
void_part=find(obst);
porosity=(size(void_part,1))/(size(obst,1)*size(obst,2));

ni=zeros(xmax,4);
ni(:,1)=(xmin:xmax)+1;
ni(:,2)=(ymin:ymax)+1;
ni(:,3)=(xmin:xmax)-1;
ni(:,4)=(ymin:ymax)-1;

% Fix near neighbours at edges (Periodic boundary condition)
ni(xmax,1)=xmin;
ni(ymax,2)=ymin;
ni(xmin,3)=xmax;
ni(ymin,4)=ymax;

%%Calculate the equilibrium distribusion
%%Identify near neighbors
ie=ni(:,1); jn=ni(:,2); iw=ni(:,3); js=ni(:,4);
Diff1(:,:)=obst(:,:)-obst(ie,:);
Diff2(:,:)=obst(:,:)-obst(:,jn);

[row1,col1]=find(Diff1==1);
x_east=transpose(row1);
y_east=transpose(col1);
[row2,col2]=find(Diff1==-1);
x_west=transpose(1+mod(row2-1+cx(2)+xmax,xmax));
y_west=transpose(1+mod(col2-1+cy(2)+ymax,ymax));

[row3,col3]=find(Diff2==1);
x_north=transpose(row3);
y_north=transpose(col3);
[row4,col4]=find(Diff2==-1);
x_south=transpose(1+mod(row4-1+cx(3)+xmax,xmax));
y_south=transpose(1+mod(col4-1+cy(3)+ymax,ymax));

obstRegion=find(obst<1);
NonObstRegion=find(obst);

%%Fluid properties (initial average density)
tauPhi=0.7; phistar=1; rhon=1;

mat=-ones(xmax,ymax);
idx=randperm(size(void_part,1),round(size(void_part,1)*sw));
fluid_rand_part=void_part(idx);
mat(fluid_rand_part)=1;

%phi0=zeros(xmax,ymax);
phi0=mat; IntWidth=4; sigma=0.003; Gamma=10; abspermeability=825;

%%Body Forces
fx=1.5*10^-5; fy=0;

%%Chemical potential parameters
alpha=0.75*sigma/(IntWidth*phistar^4);
kappa=0.5*alpha*(IntWidth*phistar)^2;

%%Wettability Spesification
angle=deg2rad(contact_angle);
beta=acos(sin(angle)^2);
coef=2*sign(pi/2-angle)*(cos(beta/3)*(1-cos(beta/3)))^0.5;
omega=coef*(2*kappa*alpha)^.5;

%%Initialize near neighbors and order parameters phi
phin=phi0(:,:);

%%Laplacian calculation for Fluids
lapPhiX=(phi0(ie,jn)+phi0(ie,js)+phi0(iw,jn)+phi0(iw,js)...
    +4*(phi0(ie,:)+phi0(iw,:)+phi0(:,jn)+phi0(:,js))-20*phi0(:,:))/6;
lapPhiY=lapPhiX;

%%Laplacians for South obstacle
x_south_east=1+mod(x_south-1+cx(2)+xmax,xmax);
x_south_west=1+mod(x_south-1+cx(4)+xmax,xmax);

y_south_north=1+mod(y_south-1+cy(3)+ymax,ymax);
y_south_north_per=1+mod(y_south_north-1+cy(3)+ymax,ymax);

lapPhiX(sub2ind(size(lapPhiX),x_south,y_south))=(phi0(sub2ind(size(phi0),x_south_east,y_south))-2*phin(sub2ind(size(phin),x_south,y_south))+phi0(sub2ind(size(phi0),x_south_west,y_south)));
lapPhiY(sub2ind(size(lapPhiY),x_south,y_south))=(6*omega/kappa+phi0(sub2ind(size(phi0),x_south,y_south_north_per))+4*phi0(sub2ind(size(phi0),x_south,y_south_north))-5*phin(sub2ind(size(phin),x_south,y_south)))/4;

%%Laplacians for North obstacle
x_north_east=1+mod(x_north-1+cx(2)+xmax,xmax);
x_north_west=1+mod(x_north-1+cx(4)+xmax,xmax);

y_north_south=1+mod(y_north-1+cy(5)+ymax,ymax);
y_north_south_per=1+mod(y_north_south-1+cy(5)+ymax,ymax);

lapPhiX(sub2ind(size(lapPhiX),x_north,y_north))=(phi0(sub2ind(size(phi0),x_north_east,y_north))-2*phin(sub2ind(size(phin),x_north,y_north))+phi0(sub2ind(size(phi0),x_north_west,y_north)));
lapPhiY(sub2ind(size(lapPhiY),x_north,y_north))=(6*omega/kappa+phi0(sub2ind(size(phi0),x_north,y_north_south_per))+4*phi0(sub2ind(size(phi0),x_north,y_north_south))-5*phin(sub2ind(size(phin),x_north,y_north)))/4;

%%Laplacians for West obstacle
x_west_east=1+mod(x_west-1+cx(2)+xmax,xmax);
x_west_east_per=1+mod(x_west_east-1+cx(2)+xmax,xmax);

y_west_south=1+mod(y_west-1+cx(5)+ymax,ymax);
y_west_north=1+mod(y_west-1+cx(3)+ymax,ymax);

lapPhiX(sub2ind(size(lapPhiX),x_west,y_west))=(6*omega/kappa+phi0(sub2ind(size(phi0),x_west_east_per,y_west))+4*phi0(sub2ind(size(phi0),x_west_east,y_west))-5*phin(sub2ind(size(phin),x_west,y_west)))/4;
lapPhiY(sub2ind(size(lapPhiY),x_west,y_west))=(phi0(sub2ind(size(phi0),x_west,y_west_south))-2*phin(sub2ind(size(phin),x_west,y_west))+phi0(sub2ind(size(phi0),x_west,y_west_north)));

%%Laplacian for the East obstacle
x_east_west=1+mod(x_east-1+cx(4)+xmax,xmax);
x_east_west_per=1+mod(x_east_west-1+cx(4)+xmax,xmax);

y_east_south=1+mod(y_east-1+cx(5)+ymax,ymax);
y_east_north=1+mod(y_east-1+cx(3)+ymax,ymax);

lapPhiX(sub2ind(size(lapPhiX),x_east,y_east))=(6*omega/kappa+phi0(sub2ind(size(phi0),x_east_west_per,y_east))+4*phi0(sub2ind(size(phi0),x_east_west,y_east))-5*phin(sub2ind(size(phin),x_east,y_east)))/4;
lapPhiY(sub2ind(size(lapPhiY),x_east,y_east))=(phi0(sub2ind(size(phi0),x_east,y_east_south))-2*phin(sub2ind(size(phin),x_east,y_east))+phi0(sub2ind(size(phi0),x_east,y_east_north)));

%%Local values of the phase and the chemical potential
muphin=4*alpha*phin.*(phin.^2-phistar^2)-kappa*(lapPhiX+lapPhiY);

%%Macroscopic Velocity values
ux=zeros(xmax,ymax);
uy=zeros(xmax,ymax);

%%Equilibrium distribution for zero initial velocity(Set initial density distribution of g PDF)
for a=1:9
    cu=cx(a).*ux+cy(a).*uy;
    if a==1
        gIn(a,:,:)=phin(:,:)-5/3.*Gamma.*muphin+3*w(a)*phin(:,:).*cu;
    else
        gIn(a,:,:)=3*w(a)*(Gamma.*muphin(:,:)+phin(:,:).*cu);
    end
end

%%Equilibrium distribution for zero initial velocity(Set initial density distribution of f PDF)
for a=1:9
    
    if a==1
        fIn(a,:,:)=(4/9).*(9/4.*rhon-15/4.*(phin.*muphin+1/3.*rhon));
    else
        fIn(a,:,:)=3*w(a)*(phin.*muphin+1/3.*rhon);
    end
end

%%Initialize pressure and velocity arrays
phi=zeros(xmax,ymax);

%% B: Main Code
 counter_frame=1;
 x=1:xmax; y=1:ymax;

 tic
for cycle=1:maxT
    phi2=phi;
    %%1: Calculation of macroscopic variables
    
    %%Order parameter
    phi=reshape(sum(gIn),xmax,ymax);
    
    %%Local values of the order parameter
    phin=phi(:,:);
    
    %%Laplacian for Fluids
    lapPhiX=(phi(ie,jn)+phi(ie,js)+phi(iw,jn)+phi(iw,js)+...
        +4*(phi(ie,:)+phi(iw,:)+phi(:,jn)+phi(:,js))-20*phin(:,:))/6;
    lapPhiY=lapPhiX;

    %%Laplacians for South obstacle
    x_south_east=1+mod(x_south-1+cx(2)+xmax,xmax);
    x_south_west=1+mod(x_south-1+cx(4)+xmax,xmax);

    y_south_north=1+mod(y_south-1+cy(3)+ymax,ymax);
    y_south_north_per=1+mod(y_south_north-1+cy(3)+ymax,ymax);

    lapPhiX(sub2ind(size(lapPhiX),x_south,y_south))=(phi(sub2ind(size(phi),x_south_east,y_south))-2*phin(sub2ind(size(phin),x_south,y_south))+phi(sub2ind(size(phi),x_south_west,y_south)));
    lapPhiY(sub2ind(size(lapPhiY),x_south,y_south))=(6*omega/kappa+phi(sub2ind(size(phi),x_south,y_south_north_per))+4*phi(sub2ind(size(phi),x_south,y_south_north))-5*phin(sub2ind(size(phin),x_south,y_south)))/4;

    %%Laplacians for North obstacle
    x_north_east=1+mod(x_north-1+cx(2)+xmax,xmax);
    x_north_west=1+mod(x_north-1+cx(4)+xmax,xmax);

    y_north_south=1+mod(y_north-1+cy(5)+ymax,ymax);
    y_north_south_per=1+mod(y_north_south-1+cy(5)+ymax,ymax);

    lapPhiX(sub2ind(size(lapPhiX),x_north,y_north))=(phi(sub2ind(size(phi),x_north_east,y_north))-2*phin(sub2ind(size(phin),x_north,y_north))+phi(sub2ind(size(phi),x_north_west,y_north)));
    lapPhiY(sub2ind(size(lapPhiY),x_north,y_north))=(6*omega/kappa+phi(sub2ind(size(phi),x_north,y_north_south_per))+4*phi(sub2ind(size(phi),x_north,y_north_south))-5*phin(sub2ind(size(phin),x_north,y_north)))/4;

    %%Laplacians for West obstacle
    x_west_east=1+mod(x_west-1+cx(2)+xmax,xmax);
    x_west_east_per=1+mod(x_west_east-1+cx(2)+xmax,xmax);

    y_west_south=1+mod(y_west-1+cx(5)+ymax,ymax);
    y_west_north=1+mod(y_west-1+cx(3)+ymax,ymax);

    lapPhiX(sub2ind(size(lapPhiX),x_west,y_west))=(6*omega/kappa+phi(sub2ind(size(phi),x_west_east_per,y_west))+4*phi(sub2ind(size(phi),x_west_east,y_west))-5*phin(sub2ind(size(phin),x_west,y_west)))/4;
    lapPhiY(sub2ind(size(lapPhiY),x_west,y_west))=(phi(sub2ind(size(phi),x_west,y_west_south))-2*phin(sub2ind(size(phin),x_west,y_west))+phi(sub2ind(size(phi),x_west,y_west_north)));

    %%Laplacian for the East obstacle
    x_east_west=1+mod(x_east-1+cx(4)+xmax,xmax);
    x_east_west_per=1+mod(x_east_west-1+cx(4)+xmax,xmax);

    y_east_south=1+mod(y_east-1+cx(5)+ymax,ymax);
    y_east_north=1+mod(y_east-1+cx(3)+ymax,ymax);

    lapPhiX(sub2ind(size(lapPhiX),x_east,y_east))=(6*omega/kappa+phi(sub2ind(size(phi),x_east_west_per,y_east))+4*phi(sub2ind(size(phi),x_east_west,y_east))-5*phin(sub2ind(size(phin),x_east,y_east)))/4;
    lapPhiY(sub2ind(size(lapPhiY),x_east,y_east))=(phi(sub2ind(size(phi),x_east,y_east_south))-2*phin(sub2ind(size(phin),x_east,y_east))+phi(sub2ind(size(phi),x_east,y_east_north)));

    %%Gradients
    gradPhiX=(4*(phi(ie,:)-phi(iw,:))+phi(ie,jn)-phi(iw,js)...
        +phi(ie,js)-phi(iw,jn))/12;
    gradPhiY=(4*(phi(:,jn)-phi(:,js))+phi(ie,jn)-phi(iw,js)...
        -phi(ie,js)+phi(iw,jn))/12;
    
    %%Gradients for wet-defined obstacle
    %%Gradient of South wall
    gradPhiX(sub2ind(size(gradPhiX),x_south,y_south))=(phi(sub2ind(size(phi),x_south_east,y_south))-phi(sub2ind(size(phi),x_south_west,y_south)))/2;
    gradPhiY(sub2ind(size(gradPhiY),x_south,y_south))=-omega/kappa;
   
    %%Gradient of North wall
    gradPhiX(sub2ind(size(gradPhiX),x_north,y_north))=(phi(sub2ind(size(phi),x_north_east,y_north))-phi(sub2ind(size(phi),x_north_west,y_north)))/2;
    gradPhiY(sub2ind(size(gradPhiY),x_north,y_north))=-omega/kappa;
    
    %%Gradient of west wall
    gradPhiX(sub2ind(size(gradPhiX),x_west,y_west))=-omega/kappa;
    gradPhiY(sub2ind(size(gradPhiY),x_west,y_west))=(phi(sub2ind(size(phi),x_west,y_west_south))-phi(sub2ind(size(phi),x_west,y_west_north)))/2;
    
    %%Gradient of East wall
    gradPhiX(sub2ind(size(gradPhiX),x_east,y_east))=-omega/kappa;
    gradPhiY(sub2ind(size(gradPhiY),x_east,y_east))=(phi(sub2ind(size(phi),x_east,y_east_south))-phi(sub2ind(size(phi),x_east,y_east_north)))/2;

    %%Chemical potential
    muPhi=4*alpha*phin.*(phin.^2-phistar^2)-kappa*(lapPhiX+lapPhiY);

    %%Hydrodynamics
    %%Local density values
    rhon=reshape(sum(fIn),xmax,ymax); 
    tauRho=3*(matdynvisc1).^((rhon+phi)./(2*rhon)).*(matdynvisc2).^((rhon-phi)./(2*rhon))+0.5;
    
    for iii=1:9
        taurho(iii,:,:)=tauRho(:,:);
    end
    
    p=alpha*(3*phi.^4-2*phistar*phi.^2-phistar^4)-kappa*(phi.*(lapPhiX+lapPhiY)+0.5*(gradPhiX.^2+gradPhiY.^2))+rhon/3;
    
    %%Interfacial force+ Body force
    sFx=muPhi.*gradPhiX+fx;
    sFy=muPhi.*gradPhiY+fy;
    
    %%Velocity field at each node
    ux=(reshape((cx*reshape(fIn,9,xmax*ymax)),xmax,ymax))./rhon;  
    uy=(reshape((cy*reshape(fIn,9,xmax*ymax)),xmax,ymax))./rhon;
        
    %%Collision for distribution function g
    for a=1:9
        newx=1+mod(x-1+cx(a)+xmax,xmax);
        newy=1+mod(y-1+cy(a)+ymax,ymax);
        cu=cx(a).*ux(:,:)+cy(a).*uy(:,:);
        if a==1
            gEq(a,:,:)=phin(:,:)-(5/3).*Gamma.*muPhi+3*w(a)*phin(:,:).*cu;
        
        else
            gEq(a,:,:)=3*w(a)*(Gamma.*muPhi(:,:)+phin(:,:).*cu);
        end
        gout(a,:,:)=gIn(a,:,:)+(gEq(a,:,:)-gIn(a,:,:))/tauPhi;
    end
    
    %%Collision for distribution function f
    for a=1:9
        
        cu=3.*(cx(a).*ux(:,:)+cy(a).*uy(:,:));
        su=w(a).*rhon(:,:).*(cu+0.5*(cu.*cu)-1.5*(ux(:,:).^2+uy(:,:).^2));
        
        if a==1
            fEq(a,:,:)=w(a).*(9/4.*rhon(:,:)-15/4.*(phin(:,:).*muPhi(:,:)+1/3.*rhon(:,:)))+su;
        else
            fEq(a,:,:)=3*w(a).*(phin(:,:).*muPhi(:,:)+1/3.*rhon(:,:))+su;
        end

        FFx(a,:,:)=cx(a).*(sFx);
        FFy(a,:,:)=cy(a).*(sFy);
        fout(a,:,:)=fIn(a,:,:)-(1./taurho(a,:,:)).*(fIn(a,:,:)-fEq(a,:,:))+3*w(a)*(FFx(a,:,:)+FFy(a,:,:));
 
    end
       
    %% C: Applying the boundary condition for fi and gi

    %%Collision of defined-obstacle
    for i=1:9
        fout(i,obstRegion)=fIn(opp(i),obstRegion);
        gout(i,obstRegion)=gIn(opp(i),obstRegion);
    end

    %%5: Streamin step
    for a=1:9
        fIn(a,:,:)=circshift(fout(a,:,:),[0,cx(a),cy(a)]);
        gIn(a,:,:)=circshift(gout(a,:,:),[0,cx(a),cy(a)]);
    end

    %% D: Visualization:
    if mod(cycle,tplot)==1
        cycle
        u = reshape(sqrt(ux.^2+uy.^2),xmax,ymax);
        phi(obstRegion) = nan;
        ux(obstRegion)=nan;
        imagesc(phi');

        %axis equal 
        colorbar;
        F(counter_frame)=getframe;
        counter_frame=counter_frame+1;
        drawnow;
        uxp1=0; lp1=0; uxn1=0; ln1=0;
        uxp2=0;lp2=0;uxn2=0;ln2=0;
        uxp3=0;lp3=0;uxn3=0;ln3=0;
        
        for yyy=ymin:ymax
            if obst(xmin,yyy)==1
                if phi(xmin,yyy)>0
                    uxp1=uxp1+ux(xmin,yyy);
                    lp1=lp1+1;
                elseif phi(xmin,yyy)<0
                    uxn1=uxn1+ux(xmin,yyy);
                    ln1=ln1+1;
                end
            end
        end
        uxp1_mean=uxp1/lp1; uxn1_mean=uxn1/ln1;
        kabs_p1=lp1/ymax*uxp1_mean*rho1*dynvisc1/fx
        kabs_n1=ln1/ymax*uxn1_mean*rho2*dynvisc2/fx
        
        for yyy=ymin:ymax
            if obst(xmax/2,yyy)==1
                if phi(xmax/2,yyy)>0
                    uxp2=uxp2+ux(xmax/2,yyy);
                    lp2=lp2+1;
                elseif phi(xmax/2,yyy)<0
                    uxn2=uxn2+ux(xmax/2,yyy);
                    ln2=ln2+1;
                end
            end
        end
        
        uxp2_mean=uxp2/lp2; uxn2_mean=uxn2/ln2;
        kabs_p2=lp2/ymax*uxp2_mean*rho1*dynvisc1/fx
        kabs_n2=ln2/ymax*uxn2_mean*rho2*dynvisc2/fx

        for yyy=ymin:ymax
            if obst(xmax,yyy)==1
                if phi(xmax,yyy)>0
                    uxp3=uxp3+ux(xmax,yyy);
                    lp3=lp3+1;
                elseif phi(xmax,yyy)<0
                    uxn3=uxn3+ux(xmax,yyy);
                    ln3=ln3+1;
                end
            end
        end
        
        uxp3_mean=uxp3/lp3; uxn3_mean=uxn3/ln3;
        kabs_p3=lp3/ymax*uxp3_mean*rho1*dynvisc1/fx;
        kabs_n3=ln3/ymax*uxn3_mean*rho2*dynvisc2/fx;
        diferent=abs(phi2(NonObstRegion)-phi(NonObstRegion));
        eps=max(diferent(:));
        
       toc
    end

    if eps <= conv
        break
    end
    
    
end





