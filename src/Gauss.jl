module Gauss
using DynamicPolynomials

export LGS, LGSvar
export elim_with!
export rückwärtsauflösen!
export gauss!
export pivot_max, pivot_first
export gerade 

abstract type AbstractLGS end

mutable struct LGS <: AbstractLGS
    A::Matrix{Float64}
    b::Vector{Float64}
    rows::Int
    cols::Int
    function LGS(A::Matrix{T}, b::Vector{S}) where {S, T}
        A = convert(Matrix{Float64}, A)
        b = convert(Vector{Float64}, b)
        n, m = size(A)
        if n < m
            Aa = zeros(m,m)
            Aa[1:n,1:m] = A
            A = Aa
            bb = zeros(m)
            bb[1:n] = b
            b = bb
            n = m
        end
        new(A, b, n, m)
    end
end

function LGS(Ab::Matrix{T}) where T
    m,n = size(Ab)
    return LGS(Ab[1:m,1:n-1], Ab[1:m, n])
end

rows(lgs::AbstractLGS) = lgs.rows
cols(lgs::AbstractLGS) = lgs.cols

function Base.show(io::IO, lgs::AbstractLGS)
    for i = 1:rows(lgs)
        for j = 1:cols(lgs)
            print(io, lgs.A[i,j])
            print(io, "\t ")
        end
        print(io, "| " )
        println(io, lgs.b[i])
    end
    return
end

function row_times!(lgs::AbstractLGS, row::Int, factor; verb = false)
    lgs.A[row,:] = factor.*copy(lgs.A[row,:])
    lgs.b[row] = factor.*copy(lgs.b[row])
    if verb
        println("Multiply row $row by $factor")
        print(lgs.A[row,:])
        print("|")
        println(lgs.b[row])
    end
    return lgs
end

function plus_row_to_row!(lgs::AbstractLGS, row1::Int, row2::Int; verb = false)
    lgs.A[row2,:] += copy(lgs.A[row1,:])
    lgs.b[row2] += copy(lgs.b[row1])
    if verb
        println("Add row $row1 to $row2")
        println(lgs)
    end
    return lgs
end

function plus_row_to_row_times!(lgs::AbstractLGS, row1::Int, row2::Int, factor::Number; verb = false)
    row_times!(lgs, row2, factor; verb = verb)
    plus_row_to_row!(lgs, row1, row2; verb = verb)
    return lgs
end

function swap_row_with!(lgs::AbstractLGS, row1::Int, row2::Int; verb = false)
    foo = copy(lgs.A[row1,:])
    bar = copy(lgs.b[row1])
    lgs.A[row1,:] = copy(lgs.A[row2,:])
    lgs.b[row1] = copy(lgs.b[row2])
    lgs.A[row2,:] = foo
    lgs.b[row2] = bar
    if verb
        println("Swap row $row1 with row $row2")
    end
    return lgs
end

function elim_with!(lgs::AbstractLGS, row1::Int, col::Int, row2::Int; verb = false)
    if lgs.A[row1,col] != 0
        plus_row_to_row_times!(lgs, row2, row1, -lgs.A[row2,col]/lgs.A[row1,col];verb = verb)
    end
    return lgs
end

abstract type PivotMode end
struct pivot_max <: PivotMode end
struct pivot_first <: PivotMode end


function pivot_element!(lgs::AbstractLGS, row::Int, col::Int; verb = false, pivot = pivot_max)

    if pivot == pivot_max
        p, p_ind = findmax(abs.(lgs.A[row:rows(lgs),col]))
    elseif pivot == pivot_first
        p_ind = findfirst(abs.(lgs.A[row:rows(lgs),col]).>0)

        if p_ind isa Nothing
            p = 0
        else
            p = lgs.A[row+p_ind-1, col]
        end
    end

    if isapprox(p, 0)
        return false
    else
    if verb
        println("Pivot for col $col is $p in row $(p_ind+row-1)")
    end
        swap_row_with!(lgs, row, p_ind+row-1; verb = verb)
        return true
    end
end

function make_zeros_below!(lgs::AbstractLGS, row::Int, col::Int; verb = false)
    for i = row+1:rows(lgs)
        elim_with!(lgs, i, col, row; verb = verb)
    end
    return lgs
end

function gauss!(lgs::AbstractLGS;  verb = false, pivot = pivot_max)
    println("")
    println("Betrachte:")
    println(lgs)
    for i = 1:rows(lgs) 
        for j = i:cols(lgs)
        if pivot_element!(lgs,i,j; verb = verb, pivot = pivot)
            make_zeros_below!(lgs,i,j; verb = verb)
            break
        end
        end
    end
    println("Gauss done:")
    println(lgs)
    return lgs
end

function rückwärtsauflösen!(lgs::AbstractLGS; verb = false)
    for i = 2:rows(lgs)
        for j = 1 : min(cols(lgs), i-1)
            @assert lgs.A[i,j] == 0
        end
    end

    for i = cols(lgs):-1:1
        row_times!(lgs, i, 1/lgs.A[i,i]; verb = verb)
        for j = i-1:-1:1
            elim_with!(lgs, j, i, i, verb = verb)
        end
    end
    println("Aufgelöst:")
    println(lgs)
    return lgs 
end

mutable struct LGSvar <: AbstractLGS
    A::Matrix{Float64}
    b::Vector{Polynomial{true,Float64}}
    rows::Int
    cols::Int
    function LGSvar(A::Matrix{T}, b::Vector{S}) where {S, T}
        A = convert(Matrix{Float64}, A)
        b = convert(Vector{Polynomial{true,Float64}}, b)
        n,m = size(A)
        if n < m
            Aa = zeros(m,m)
            Aa[1:n,1:m] = A
            A = Aa
            bb = zeros(m)
            bb[1:n] = b
            b = bb
            n = m
        end
        @assert n == length(b)
        new(A, b, n, m)
    end
end

LGS(A::Matrix{Float64}, b::Vector{Polynomial{true, Float64}}) = LGSvar(A,b)

function gerade(lgs::LGS; var = "t")
    for k = cols(lgs):rows(lgs)
        @assert all(isapprox.(lgs.A[k,:], 0))
    end

    eval(Meta.parse("foo = first(@polyvar "*var*")"))
    A = copy(lgs.A)
    b = convert.(Polynomial{true, Float64},copy(lgs.b))
  
    if cols(lgs)<rows(lgs)
        A = [A;zeros(cols(lgs)-rows(lgs),cols(lgs))]
        b = [b;zeros(cols(lgs)-rows(lgs))]
    end
    A[cols(lgs),rows(lgs)] = 1
    b[cols(lgs)] = foo
    

    return rückwärtsauflösen!(LGSvar(A,b))
end

end # module Gauss
