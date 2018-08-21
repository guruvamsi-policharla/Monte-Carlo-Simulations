#Includes
include("ising_aux.jl")
using Distributions
using Plots

#Code
Tmin = 0.01
Tchange = 0.1
Tmax = 5.0
mcs = 1000
M = 4
N = 4


lat = initialise(M,N)

norm=(1.0/float(M*N))

Temperature = Tmin:Tchange:Tmax

E_vec = zeros(length(Temperature),2)
M_vec = zeros(length(Temperature),2)
acc_vec = zeros(length(Temperature),2)

E_jack = zeros(mcs,1)
M_jack = zeros(mcs,1)
acc_jack = zeros(mcs,1)


count = 1
for T in Temperature

    transient_results(lat,1000,T)
    Mag = total_mag(lat)
    E = total_energy(lat)

    etot=0;
    mabstot=0;
    acc_rat = 0;

    for i in 1:mcs
        for j in 1:M*N
            x = rand(1:M)
            y = rand(1:N)
            E_0 = energy_pos(x,y,lat)
            Mag_0 = lat[x,y]
            if(test_flip(x,y,lat,T))
                acc_rat = acc_rat + 1
                E = E + energy_pos(x,y,lat) - E_0
                Mag = Mag + lat[x,y] - Mag_0
            end
        end
        mabs = sqrt(vecdot(Mag,Mag))

        E_jack[i] = E
        M_jack[i] = mabs
        acc_jack[i] = acc_rat

        etot= etot + E
        mabstot = mabstot + mabs
    end
    E_jack = E_jack*norm
    M_jack = M_jack*norm
    acc_jack = acc_jack*norm

    M_vec[count,1], M_vec[count,2] = jackknife(M_jack)
    E_vec[count,1], E_vec[count,2] = jackknife(E_jack)
    acc_vec[count,1], acc_vec[count,2] = jackknife(acc_jack)
    println(T," ",E_vec[count,1])
    count = count + 1
end

#jackknife analysis
