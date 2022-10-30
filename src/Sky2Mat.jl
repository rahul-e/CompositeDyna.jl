function Sky2Mat(v::Vector{Float64}, nbig::Int64, nsky::Int64, cht::Vector{Int64})

    G=fill(0.0,(nbig,nbig))
    k=1
    for i=1:nbig
        for j=i:-1:1
            ii=i-cht[i]
            if(k<=nsky && j>=ii)
                G[j,i]=v[k]
                k=k+1
            end
        end
    end
    return G  
end