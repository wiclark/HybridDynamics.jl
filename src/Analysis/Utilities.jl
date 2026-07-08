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

# Isolate only the jump transitions
function extract_jumps(sol::AbstractHybridSolution)
    states = sol.x

    T = typeof(states[1])
    data_list = T[]
    t_list = Float64[]

    #Create NaN container to break the plot lines between different jumps
    nan_state = fill!(similar(states[1]), NaN)

    for i in 1:(length(states)-1)
        #A jump happens when times doesnt advance
        if sol.t[i] == sol.t[i+1]
            #Add pre jump state
            push!(data_list, states[i])
            push!(t_list, sol.t[i])

            #Add post jump state
            push!(data_list, states[i+1])
            push!(t_list, sol.t[i+1])

            #Insert NaN break so it doesnt connect to the next jump event
            push!(data_list, nan_state)
            push!(t_list, NaN)
        end
    end
    return t_list, data_list
end