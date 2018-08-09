#Includes
include("ising_aux.jl")
using Distributions
using Plots

#Code
Tmin = 0.001
Tchange = 0.01
Tmax = 5.0
mcs = 100000
N = 4
acc_rat = 0;

lat = initialise(N)

norm=(1.0/float(mcs*N*N))

Temperature = Tmin:Tchange:Tmax
E_vec = zeros(length(Temperature),1)
Mabs_vec = zeros(length(Temperature),1)
acc_vec = zeros(length(Temperature),1)
count = 1
for T in Temperature
    transient_results(lat,1000,T)
    M = total_mag(lat)
    Mabs = abs.(M)
    E = total_energy(lat)

    etot=0;
    mabstot=0;

    for i in 1:mcs
        for j in 1:N*N
            x = rand(1:N)
            y = rand(1:N)
            E_0 = energy_pos(x,y,lat)
            if(test_flip(x,y,lat,T))
                acc_rat = acc_rat + 1
                E = E + energy_pos(x,y,lat) - E_0
            end
        end
        etot= etot + E
    end
    acc_rat = acc_rat*norm
    E_avg=etot*norm;
    E_vec[count] = E_avg
    acc_vec[count] = acc_rat
    count = count + 1
    println(T," ",E_avg)
end

#jackknife analysis
