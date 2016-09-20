#Greens function 
Gn<-function(x1, t1, mu){
   # global D muN muP gamma  
    (1/sqrt(4*pi*D*(t1)))*exp(-mu*(t1)-(x1-gamma*(t1))^2/(4*D*(t1)));
    }

#forcing function
fset<-function(t0,r1,r2,tau){
   # global f0 
    if (tau <= t0){forcing = f0*exp(r1*(tau-t0))}else{forcing = f0*exp(r2*(tau-t0))}
    forcing}
   