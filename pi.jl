@everywhere function compute_pi(N::Int)
    """@everywhere tells all cores about the existence of compute_pi(N)"""
    count = 0;
    for i = 1:N
        x = rand()*2 - 1;
        y = rand()*2 - 1;
        r2 = x.^2 + y.^2;
        if r2< 1.0
            count = count + 1;
        end
    end
    return count/N*4.0;
end

function parallel_pi_computation(N::Int; ncores::Int=8)
    """
    Compute pi in parallel, over ncores cores, with a Monte Carlo simulation throwing N total darts
    """

    # compute sum of pi's estimated among all cores in parallel
    sum_of_pis = @parallel (+) for i=1:ncores
        compute_pi(ceil(Int, N / ncores))
    end

    return sum_of_pis / ncores  # average value
end
