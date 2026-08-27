#Plotting Lines Expelled!
"""
    split_jumps(sol::AbstractHybridSolution)

Inserts NaNs to separate jumps in plots
"""
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
"""
    extract_jumps(sol::AbstractHybridSolution)

Isolate the times/states when a jump occurs.
"""
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

# Input time and get out how many jumps have occurred up to that time
"""
    jump_count(sol::AbstractHybridSolution, t::Real)

Counts how many jumps/events have occured up to time t.
"""
function jump_count(sol::AbstractHybridSolution, t::Real)
    # Ensure the requested time is within the simulation
    if isempty(sol.event_times) || t < first(sol.t)
        return 0
    end

    # searchsortedlast returns the index of the last element in the array
    return searchsortedlast(sol.event_times, t)
end

# Returns the interval '(t_start, t_end)' that contains the jump k
"""
    jump_interval(sol::AbstractHybridSolution, k::Int)

Returns the time interval between events k and k+1.
"""
function jump_interval(sol::AbstractHybridSolution, k::Int)
    if k<0
        throw(ArgumentError("Jump count `k` must be non-negative."))
    end

    n_jumps = length(sol.event_times)

    if k > n_jumps
        throw(ArgumentError("Jump count `k` exceeds the number of jumps in the solution."))
    end

    # if k = 0, the intervals starts at the sims initial time
    t_start = k == 0 ? first(sol.t) : sol.event_times[k]

    # if k is max number of jumps reached, the interval ends at the sim final time. 
    t_end = k == n_jumps ? last(sol.t) : sol.event_times[k+1]

    return (t_start, t_end)
end