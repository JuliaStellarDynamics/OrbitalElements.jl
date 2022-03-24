# what is the optimal spacing to build an interpolation table?
# definitely log radius, with some base,
# but what about eccentricity? does it help to pile near the edges?


function make_frequency_grid(potential::Function,dpotential::Function,ddpotential::Function,zerofreq::Float64,
                             amin::Float64,amax::Float64,res::Int64)
    avals = LinRange(amin,amax, res)
    evals = LinRange(0., 1., res)
    agrid  = zeros(Float64, res,res)
    egrid  = zeros(Float64, res,res)
    O1grid = zeros(Float64, res,res)
    O2grid = zeros(Float64, res,res)
    alphagrid = zeros(Float64, res,res)
    betagrid = zeros(Float64, res,res)
    dO1dagrid = zeros(Float64, res,res)
    dO2dagrid = zeros(Float64, res,res)
    dO1degrid = zeros(Float64, res,res)
    dO2degrid = zeros(Float64, res,res)

    for a=1:res
        semia = 10^avals[a]

        for b=1:res
            ecc = evals[b]
            rp,ra = rpra_from_ae(semia,ecc)
            O1grid[b,a],O2grid[b,a] = compute_frequencies_henon_ae(potential,dpotential,ddpotential,semia,ecc)

            #O1grid[b,a],O2grid[b,a] = compute_frequencies_henon_ae(isochrone_psi,isochrone_dpsi_dr,isochrone_ddpsi_ddr,semia,ecc)
            #O1grid[b,a],O2grid[b,a] = compute_frequencies_henon_ae(plummer_psi,plummer_dpsi_dr,plummer_ddpsi_ddr,semia,ecc)
            alphagrid[b,a] = O1grid[b,a]/zerofreq
            #alphagrid[b,a] = O1grid[b,a]/plummer_Omega0()
            betagrid[b,a] = O2grid[b,a]/O1grid[b,a]

            agrid[b,a] = semia
            egrid[b,a] = ecc
        end
    end


    # create empirical diff grids
    for a=1:res    # semi-a values
        for b=1:res # eccentricity values
            if a==res
                dO1dagrid[b,a] = (O1grid[b,a]-O1grid[b,a-1])/(agrid[b,a]-agrid[b,a-1])
                dO2dagrid[b,a] = (O2grid[b,a]-O2grid[b,a-1])/(agrid[b,a]-agrid[b,a-1])
            else
                dO1dagrid[b,a] = (O1grid[b,a+1]-O1grid[b,a])/(agrid[b,a+1]-agrid[b,a])
                dO2dagrid[b,a] = (O2grid[b,a+1]-O2grid[b,a])/(agrid[b,a+1]-agrid[b,a])
            end

            if b==res
                dO1degrid[b,a] = (O1grid[b,a]-O1grid[b-1,a])/(egrid[b,a]-egrid[b-1,a])
                dO2degrid[b,a] = (O2grid[b,a]-O2grid[b-1,a])/(egrid[b,a]-egrid[b-1,a])
            else
                dO1degrid[b,a] = (O1grid[b+1,a]-O1grid[b,a])/(egrid[b+1,a]-egrid[b,a])
                dO2degrid[b,a] = (O2grid[b+1,a]-O2grid[b,a])/(egrid[b+1,a]-egrid[b,a])
            end
        end
    end

    return avals,evals,agrid,egrid,O1grid,O2grid,dO1dagrid,dO1degrid,dO2dagrid,dO2degrid
end
