function ColumnH(cht::Int64, nd::Vector{Int64}, ned::Int64, neq::Int64)
    # cht - column height refers to the number of elements 
    # in each column below the skyline and above the leading 
    # diagonal excluding the diagonal element of the global stiffness matrix
    # 
    # nd - nodal ordinate (y-coordinate)
    # ned - degree of freedom in an element
    # neq - number of equations

    # C. S. Krishnamoorthy Pg. 184
    # Set ls equal to the smallest equation number of all the
    # d.o.f of the element
    # 
    ls::Int64=15000
    tmp::Int64=0

    for i=1:ned
        if (nd[i]!=0) 
            if((nd[i]-ls)>=0)
                ls=nd[i]
            end
        end
    end

    for i=1:ned
        tmp=nd[i]
        if(tmp!=0)
            me=tmp-ls
            if(me > cht[tmp])
                cht[tmp]=me
            end
        end
    end
    
end

function CadNum(cht::Int64, nds::Vector{Int64}, neq::Int64, nsky::Int64,
    mband::Int64)
    # cht - column height refers to the number of elements
    # nds - address of the diagonal elements
    # neq - number of equations
    # neq1 = neq+1
    # Compute address of diagonal elements
    # when column height is known C. S. Krishnamoorthy Pg. 186

    for i=1:(neq+1)
        nds[i]=0
    end
    nds[1]=1
    nds[2]=2
    mband=0
    if(neq!=1)
        for i=2:neq
            if(cht[i]>mband)
                mband=cht[i]
            end
            nds[i+1]=nds[i]+cht[i]+1
        end
    end

    mband=mband+1
    nsky=nds[neq+1]-1
    
    return
end

function Passem(sk::Vector{Float64}, ek::Vector{Float64}, nds::Vector{Int64}, nd::Vector{Int64}, ned::Int64, neq1::Int64, nsky::Int64, nued::Int64)
    sk::Vector{Float64}(undef, nsky)
    nds::Vector{Int64}(undef, neq1)
    nd::Vector{Int64}(undef, ned)
    ek::Array{Float64}(undef, nued, nued)
    ii::Int64
    jj::Int64
    mi::Int64
    # Passem function from Krishnamoorthy Pg. 600
    for i=1:ned
        ii=nd[i]

        if(ii!=0)
            for j=1:ned
                jj=nd[j]

                if(jj!=0)
                    mi=nds[jj]
                    ij=jj-ii

                    if(ij>=0)
                        kk=mi+jj
                        sk[kk]=sk[kk]+ek[i,j]
                    end
                end
            end
        end
    end
    return sk
end


            



    

