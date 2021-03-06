include("ising_aux.jl")
using Plots
using Distributions
Tmin = 0.001
Tchange = 0.1
Tmax = 5.0
mcs = 100000
N = 4
lat = initialise(N)

norm=(1.0/float(mcs*N*N))

Temperature = Tmin:Tchange:Tmax
E_vec = zeros(length(Temperature),1)
Mabs_vec = zeros(length(Temperature),1)
count = 1
for T in Temperature

    transient_results(lat,1000,T)
    M = total_mag(lat)
    E = total_energy(lat)

    etot=0;
    mabstot=0;

    for i in 1:mcs
        for j in 1:N*N
            x = rand(1:N)
            y = rand(1:N)
            if(test_flip(x,y,lat,T))
                lat[x,y] = -lat[x,y]
                E = E + 2*energy_pos(x,y,lat)
                M = M + 2*lat[x,y]
                Mabs = abs(M)
            end
        end
        etot= etot + E
        mabstot= mabstot + abs(M);
    end

    E_avg=etot*norm;
    Mabs_avg=mabstot*norm;
    E_vec[count] = E_avg
    Mabs_vec[count] = Mabs_avg
    count = count + 1
    println(T," ",E_avg)
end

#jackknife analysis
