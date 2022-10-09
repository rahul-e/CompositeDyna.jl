function ColumnH(cht::Int64, nd::Float64, ned::Int64, neq::Int64)
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
            if((nd[i]-ls)!=0)
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

function CadNum(cht::Int64, nds::Int64, neq::Int64, nsky::Int64,
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
        end

    else
        mband=mband+1

    nsky=nds[neq+1]-1
    return
end
    

