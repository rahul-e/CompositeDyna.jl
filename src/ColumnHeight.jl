function ColumnH(cht::Vector{Int64} , nd::Vector{Int64}, ned::Int64, neq::Int64)
		
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
    ls::Int64=150000
    tmp::Int64=0
    #cht=fill(0, neq)
    #nd=Vector{Int64}(undef, ned)
    #println("neq", neq)
    
    for i=1:ned
        if (nd[i]!=0) 
            if((nd[i]-ls)<0)
                ls=nd[i]
                #println(ls)
            end
        end
    end
    # Second loop to Calculate the column Height
    for i=1:ned
        tmp=nd[i]
        #println(i," ",nd[i])
        if(tmp!=0)
            me=tmp-ls
            if(me > cht[tmp]) # Change to greater than equalto
                cht[tmp]=me
            end
        end
    end

    return cht
    
end 

function CadNum(cht::Vector{Int64}, neq::Int64)
		
    # cht - column height refers to the number of elements in each column of the global stiffness matrix below the skyline and above the leading diagonal element.
    # nds - address of the diagonal elements
    # neq - number of equations
    # neq1 = neq+1
    # Compute address of diagonal elements
    # when column height is known C. S. Krishnamoorthy Pg. 186
    #nsky::Int64=0
    # The number of diagonal elements of the global stiffness matrix ar stored in "nds" array
    nds = fill(0, neq+1)
    
    nds[1]=1
    nds[2]=2
    mband::Int64=0
    
    if neq == 1
        mband=mband+1
        nsky = nds[neq1]-1
        return nds, nsky
    end

    for i=2:neq
            if(cht[i]>mband)
                mband=cht[i]
            end
            nds[i+1]=nds[i]+cht[i]+1
    end
    mband=mband+1
    
    nsky=nds[neq+1]-1
    return nds, nsky
end

function passem!(sk::Vector{Float64},ek::Array{Float64, 2}, nds::Vector{Int64},
    nd::Vector{Int64}, ned::Int64, neq1::Int64, nsky::Int64, nued::Int64, 	 
    )
 # Passem assembles the element stiffness matrix to global stiffness matrix in the form an 1D array SK   
    i::Int64=0
    ii::Int64=0
    j::Int64=0
    jj::Int64=0
    mi::Int64=0
    ij::Int64=0
    kk::Int64=0
    #println(nsky)
    #sk = fill(0.0, nsky)
    #println(ned)

    # PASSEM function FEA - THeory and Programming, C. S. Krishnamoorthy, Pg. 190 
    for i=1:ned
        # Set ii = eqn. number corresponding to DOF(I)
        ii=nd[i]
        # Is the DOF inactive
        if (ii == 0)
            continue
        end

        for j in 1:ned
                # Set jj = eqn. number corresponding to DOF(J)
                jj=nd[j]
                # Is the DOF inactive
                if (jj == 0)
                    continue
                end
                # Set mi as the diagonal element entry number of the 
                # jj^th column of stiffness matrix

                mi=nds[jj]
                ij=ii-jj
                    
                if (ij < 0)
                        kk=mi+jj-ii
                        sk[kk]=sk[kk]+ek[i,j]
                elseif (ij == 0)
                        kk=mi
                        sk[kk]=sk[kk]+ek[i,j]
                end
                
        end
    end

    return sk
end


            



    

