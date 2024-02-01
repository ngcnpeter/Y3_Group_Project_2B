function [output]=groupproject_IRFsimulate(amplitudes,lifetimes,acquisitiontime,irfwidth)
    bg=10; %number of background counts per second, keep at 10
    if irfwidth==0
        irfwidth=0.00000001;
    end
    
    %check that each amplitude has a corresponding lifetime
    if length(amplitudes)~=length(lifetimes)
        return
    end
    
    %create empty vector to store decay data
    puredecay=zeros(381,1);
    
    %normalise amplitudes, just in case they didn't initially sum to 1
    amplitudes=amplitudes/sum(amplitudes);
    
    %generate a multiexponential decay starting at 1 at t=0
    %using the supplied amplitudes and lifetimes
    for i=1:381
        t=(1/19)*(i-1); %each bin is (1/19) ns, starting at t=0
        for j=1:length(amplitudes)
            puredecay(i,1)=puredecay(i,1)+amplitudes(j)*exp(-t/lifetimes(j));
        end
    end
    
    %generate the IRF, centred at b
    b=10/19;
    w=irfwidth;
    for i=1:381
        t=(i-20)*(1/19);
    	irfraw(i,1)=exp(-4*log(2)*(t-b)*(t-b)/(w*w));
    end
    
    %convolute the IRF and decay and trim to 381 bins
    Iconvol=conv(puredecay,irfraw);
    Iconvol=Iconvol(1:381);
    
    %we do our measurements at 2500 counts per second
    %calculate how many fluorescence counts per second this corresponds to
    %i.e. subtract background from total counts
    fluorate=2500-bg;
    
    %calculate total number of fluorescence photons counted in measurement
    %%
    % 
    % $$e^{\pi i} + 1 = 0$$
    % 
    totalfluorescence=fluorate*acquisitiontime; 
    
    %now scale the multiexponential decay so it contains this many counts
    noiseless=totalfluorescence*Iconvol/sum(Iconvol);

    %and add on 'bg' counts per second spread evenly across all bins
    noiseless=noiseless+(bg*acquisitiontime/381);
    
    %finally add poisson noise to each bin
    noisydecay=poissrnd(noiseless);
    
    %and tidy up output with a time axis
    output=zeros(381,2);
    for i=1:381
        output(i,1)=(i-1)*(1/19);
        output(i,2)=noisydecay(i,1);
    end
end