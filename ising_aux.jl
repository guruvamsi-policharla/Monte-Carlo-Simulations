function sample_gauss(v)
    #Input is the central vector around which we flip
    var = 0.1
    v_new = [rand(Normal(v[1],var)),rand(Normal(v[2],var)),rand(Normal(v[3],var))]
    v_new = v_new/sqrt(vecdot(v_new,v_new))
    return convert(Vector,v_new)
end

function sample_uni()
    x = [rand(Uniform(-1,1)), rand(Uniform(-1,1)), rand(Uniform(-1,1))]
    x = x/sqrt(vecdot(x,x))
    return convert(Vector,x)
end

function initialise(M::Int, N::Int)
    """ Initialising the lattice with random values """
    #lat = Array{Float64, 3}(N, N)
    lat = Array{Vector{Float64},2}(M, N);
    for i = 1:M
        for j = 1:N
            lat[i,j] = sample_uni()
        end
    end
    return lat
end

function energy_pos(x, y, lat, a = [0,0,0])
    M = size(lat,1)
    N = size(lat,2)


    if(M==1 && N==2)
        if(a == [0,0,0])
            energy = -1*vecdot(lat[1,1],lat[1,2])
            return energy
        else
            energy = -1*vecdot(a,lat[x,mod(y,2)+1])
            return energy
        end
    end

    if(M!=1 && N!=1)
        if(y==N)
            up=1;
        else
            up=y+1;
        end

        if(y==1)
            down=N;
        else
            down=y-1;
        end

        if(x==1)
            left=M;
        else
            left=x-1;
        end

        if(x==M)
            right=1;
        else
            right=x+1;
        end

        if(a == [0,0,0])
            energy = -1*vecdot(lat[x,y],(lat[left,y]+lat[right,y]+lat[x,up]+lat[x,down]));
            return energy
        else
            energy = -1*vecdot(a,(lat[left,y]+lat[right,y]+lat[x,up]+lat[x,down]));
            return energy
        end
    end

    return 0
end

function test_flip(x, y, lat, T)
""" Checks whether energy allows for a flip or not """
    #a = sample_uni()
    a = sample_gauss(lat[x,y])
    de = -energy_pos(x,y,lat) + energy_pos(x,y,lat,a);

    if(de<0)
        lat[x,y] = a
        return true
    elseif(rand()<exp(-de/T))
        lat[x,y] = a
        return true
    else
        return false
    end
end

function transient_results(lat, transient::Int, T)
    """Takes lat as input and removes initial transients by running for transient number of steps"""
    M = size(lat,1)
    N = size(lat,2)
    for i = 1:transient
        for j = 1:M*N
                x = rand(1:M)
                y = rand(1:N)
                test_flip(x,y,lat,T)
        end
    end
end

function total_mag(lat)
    #returns total magnetisation vector
    return sum(lat)
end

function total_energy(lat)
    e = 0.0
    for i = 1:size(lat,1)
        for j = 1:size(lat,2)
            e = e + energy_pos(i,j,lat)
        end
    end
    return e/2
end

function exact_heisen_energy(T)
    #(e^(1/T) (-1 + T) T - e^(-1/T) T (1 + T))/(2 T sinh(1/T))
    return (exp(1/T) * (-1 + T) * T - exp(-1/T) * T * (1 + T))/(2*T*sinh(1/T))
end

function jackknife(vec)
    s = sum(vec)
    n = length(vec)
    vec_jack = (s - vec)/(n-1)
    jack_avg = sum(vec_jack) / n

    jack_err = sqrt(sum((vec_jack-jack_avg).^2) * (n-1)/n)
    #jack_err = sqrt((n-1)*(jack_err - jack_avg.^2))
    return jack_avg,jack_err
end
