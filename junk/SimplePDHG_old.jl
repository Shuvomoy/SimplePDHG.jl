module SimplePDHG

include("./Types.jl")

export LP_Data, PDHG_settings, PDHG_state, PDHG_iteration!

include("./Utils.jl")

export project_nonnegative!, project_nonnegative, tolerance_LP

## The solver function

function PDHG_solver(problem::LP_Data, setting::PDHG_settings)
    
    # this is the function that the end user will use to solve a particular problem, internally it is using the previously defined types and functions to run 
    # PDHG algorithm
    # create the intial state
    state = PDHG_state(problem)

    tolerance_current =  tolerance_LP(problem.A, problem.b, problem.c, state.x, state.y)
    
    ## time to run the loop
    while  (state.k < setting.maxit) &&  tolerance_current > setting.tol
        # compute a new state
        state =  PDHG_iteration!(problem, state)
        tolerance_current =  tolerance_LP(problem.A, problem.b, problem.c, state.x, state.y)
        # print information if verbose = true
        if setting.verbose == true
            if mod(state.k, setting.freq) == 0
                @info "iteration = $(state.k) | obj val = $(problem.c'*state.x) | optimality_tolerance = $(tolerance_current)"
            end
        end
    end
    
    # print information regarding the final state
    
    @info "iteration = $(state.k) | obj val = $(problem.c'*state.x) | optimality_tolerance = $(tolerance_current)"

    return state
    
end

export PDHG_solver

end # module
