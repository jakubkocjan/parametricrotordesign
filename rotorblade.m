
clc
clearvars

R=7.85                  % rotor radius [m]
c=0.44                  % rotor chord [m]
b=4                     % number of blades [m]
rpm=255                 % rotor angular velocity [RPM]
a=6.47                  % lift curve slope [rad]
ro=1.225                % air density [kg/m^3]
G=6400 % weight
vc=0
om=rpm*pi/30  

collective=0

x=[0.21 0.275 0.36 0.45 0.54 0.63 0.72 0.81 0.9 0.95 0.975 1]  % blade coordinates
relc=[0.732 0.943 1 1 1 1 1.0955 1.05 1.05 0.827 0.664 0.5]    % chord dimensions
twist=[5.34 5.15 4.12 3.03	1.94 0.85 -0.24	-1.33 -2.42 -3.03 -3.33 -3.63] % twist value

 r=@(x) R*x                             % radius funciotn

p=polyfit(x,relc,6)
p2=polyfit(x,twist+collective,2)

chord=@(x) polyval(p,x)*c                % chord function
pitch=@(x) polyval(p2,x)                 % twist function

om=rpm*pi/30                             % rotor angular velocit [rad/s]
sig=@(x) (b.*chord(x))/(pi.*R)        % rotor soilidity

vi= (G*9.81./(2*ro*pi*R.^2)).^(1/2)
delt=@(x) atan(vi./(om.*r(x)))
collective=0
 exitFLAG = 0

B=1-chord(1)/(2*R)  

while ~exitFLAG

T=@(x) 0.5.*ro.*b*a.*om^2.*R.^2*chord(x).*r(x).*(deg2rad((pitch(x)+collective)-delt(x))).*x.^2 % thrust function

Tint=integral(T,0.21,B)            % thrust integral
Blade=Tint/4

if abs(G*9.81-Tint) < 50         % find collective loop
    exitFLAG=1
else
    collective=collective+0.01
    exitFLAG=0
end
end

Tint                             % results
Blade=Tint/4                     % results
collective                       % results
Cl=@(x) Blade./(0.5.*ro.*om^2.*R.^2*chord(x).*r(x))
Clint=integral(Cl,0.21,B)
vi
 lami=@(x) (sig(x).*a./16).*(((1+(32./(sig(x).*a)).*deg2rad(pitch(x)+collective)).^1/2)-1) % inflow ratio for hover


alf=polyfitn((rad2deg(lami(x)./x)),x,3) % inflow angle for ANSYS FLUENT
alffun=@(x) (((collective)-(polyval(alf.Coefficients,x))))  % inflow angle for ANSYS FLUENT

fplot(@(x) alffun(x),[0 1])

%influence of paramter change

%while b<6

b=3
k=0
initialtwist=5
twist=0
dtw=11
m=0
step=0.02
 k=k+1

twistang=0
twstep=-1
    while twistang>=-15
    m=m+1
    nn=0
    n=0
    c=0.3
       while nn<=dtw
    
        
           twist(1,nn+1)=initialtwist+(twistang/dtw)*nn
   
           nn=nn+1

       end

        p2=polyfit(x,twist+collective,2)
        
        pitch=@(x) polyval(p2,x)  
            while c<0.6+step
        
            n=n+1   
            chord=@(x) polyval(p,x)*c                % chord function
                           % twist function
            T=@(x) 0.5.*ro.*b*a.*om^2.*R.^2*chord(x).*r(x).*(deg2rad((pitch(x)+collective)-delt(x))).*x.^2 % thrust function
        
            Tint=integral(T,0.21,B)            % thrust integral
            Thrust(m,n)=Tint
            ch(1,n)=c
            c=c+step
            end
            tn(m,1)=twistang
            twistang=twistang+twstep
    end       
%bn(1,k)=b
%b=b+1
%end

surf(tn,ch,Thrust)

xlabel('Twist [deg]','fontsize',24), ylabel('Chord [m]','fontsize',24), zlabel('Thrust [N]','fontsize',24)
title('Main rotor thrust','fontsize',24)
set(gca, 'FontSize', 18);
view(60,40)