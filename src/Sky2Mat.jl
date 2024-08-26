function Sky2Mat(v::Vector{Float64}, nbig::Int64, nsky::Int64, cht::Vector{Int64})

    G_=fill(0.0,(nbig,nbig))
    k=1
    for col=1:nbig
        for row=col:-1:1
            limit=col-cht[col]
            if(row>=limit)
                G_[row,col]=v[k]
                k=k+1
            end
        end
    end
    G=G_+G_'- diagm(diag(G_))  

    #println("Value of k = ", k)
    println("Check for NaN in Global array ", isnan.(sum(G)))
    println("Check if all values are finite in Global array ",isfinite.(sum(G)))
    return G
end