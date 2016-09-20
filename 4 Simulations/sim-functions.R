#Greens function 
Gn<-function(x1, t1, mu){
   # global D muN muP gamma  
    (1/sqrt(4*pi*D*(t1)))*exp(-mu*(t1)-(x1-gamma*(t1))^2/(4*D*(t1)));
    }

#forcing function
fset<-function(tT, fT, r1, r2, tau){
   	Efficacy<-90
   	if (tau <= tT){
   		forcing = fT*exp(r1*(tau-tT))
    
    } else if (tau > tT & tau <= tT+Efficacy){
    	forcing = fT*exp(r2*(tau-tT))
    	
    } else if (tau > tT+Efficacy){
    		forcing = fT*exp(r2*Efficacy)*exp(r1*(tau-(tT+Efficacy)))
    	}
  
  forcing}
   