module Minc2

greet() = print("Hello World!")

export volume_plus

"""
    volume_plus(x)

Dummy volume function
# Arguments
* `x`: the input value

# Notes
* I will add notes

# Examples
```julia
julia> two = volume_plus(1)
2
```
"""
volume_plus(x) = return x+1

end # module
