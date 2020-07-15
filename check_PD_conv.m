function [primal,finish,gap_min] = check_PD_conv(primal,primal_new,dual,gap_min,i,IRGN)
gap=abs(primal_new-dual);
finish=0;
if i==1
    gap_min = gap;
end
if abs(primal-primal_new)<IRGN.tol
    % print("Terminated at iteration %d because the energy decrease in the primal problem was less than %.3e"%(i,np.abs(primal-primal_new)/(self.irgn_par["lambd"]*self.NSlice)))
    
    finish=1;
    return
end
if gap > gap_min*IRGN.stag && i>2
    
    finish=1;
    return
end

if abs(gap - gap_min)<IRGN.tol && i>2
    
    finish=1;
    return
end
primal = primal_new;
gap_min = min(gap,gap_min);




end

