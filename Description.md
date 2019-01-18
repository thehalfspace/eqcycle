Description of the slowest parts of code:
==========================================

- The main function (src/main.jl) has one time loop.

- Preconditioned Conjugate Gradient: src/PCG.jl. This is the slowest step. 
    while time < tmax:
        function PCG (src/main.jl: line 183)
            # compute initial residuals.
            
            for n = 1:4000  # number of steps for convergence
                # Element wise assembly
                function element_computations (src/PCG.jl: line 75)
                    for e = 1: no_of_elements_in_mesh
                        # compute gradients and assemble
                        local += local
                    endfor
                endfunction
            endfor
        endfunction
    endwhile

The above function has multiple nested for loops. The outer most while time loop cannot be changed. The number of steps for convergence can be improved. The element computation loop can be distributed to speed up calculations.

The element computations loops are in src/PCG.jl: line 75, line 107.


- Some other for loops which can be run in parallel, i.e., the output of one iteration does not depend on the output of previous iterations:
    src/otherFunctions.jl: line 122
    src/NRsearch.jl: line 10


