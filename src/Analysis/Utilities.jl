#Plotting Lines Expelled!
function split_jumps(sol::AbstractHybridSolution) 
    states = sol.x
    
    # Preallocate an array of exactly whatever type the state is
    T = typeof(states[1])
    data_list = T[]
    t_list = Float64[]
    
    # Create a NaN container of the exact same size and shape as the state
    nan_state = fill!(similar(states[1]), NaN)
    
    for i in 1:length(states)
        push!(data_list, states[i])
        push!(t_list, sol.t[i])
        
        # Insert the NaN break if a jump occurred (zero-time transition)
        if i < length(sol.t) && sol.t[i] == sol.t[i+1]
            push!(data_list, nan_state)
            push!(t_list, NaN)
        end
    end
    
    return t_list, data_list
end
