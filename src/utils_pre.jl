"""
Works just like `Threads.@threads`, but runs serially if `Threads.nthreads()==1` to avoid task creation overhead.
"""
macro parallel(ex)
    pex = quote
        let
            if $Threads.nthreads() == 1
                $(ex)
            else
                $Threads.@threads $(ex)
            end
        end
    end
    return esc(pex)
end 
