# Generate constitutive matrix for an orthotropic composite laminate
# with NLAYERS as number of plies
#
function composit(nlayers::Int64, ori::Vector{Float64}, zpos::Vector{Float64}, e11::Float64, 
    e22::Float64, pr12::Float64, g12::Float64, g13::Float64, g23::Float64)

    println("Laminate with number of plies = ", nlayers) #, ori, zpos, e11, e22, pr12, g12, g13, g23)

    a = Array{Float64, 2}(undef, 3, 3)
    b = Array{Float64, 2}(undef, 3, 3)
    d = Array{Float64, 2}(undef, 3, 3)
    t = Array{Float64, 2}(undef, 3, 3) #Transformation matrix
    c = Array{Float64, 2}(undef, 5, 5)
    qbar = Array{Float64, 2}(undef, 5, 5)
    qt = Array{Float64, 2}(undef, 5, 5)
    ss = Array{Float64, 2}(undef, 2, 2)
    
    # compliance matrix - this will be returned on function call
    cc = Array{Float64, 2}(undef, 8, 8)

    pr21=e22*pr12/e11


    c[1,1] = e11/(1.0-pr12*pr21)
    c[1,2] = c[1,1]*pr21 
    c[2,2] = e22/(1.0-pr12*pr21)
    c[2,1] = c[2,2]*pr12
    c[3,3] = g12
    c[4,4] = g13
    c[5,5] = g23

    #println(c)

    for l in 1:nlayers
        α=ori[l]*π/180.0
        t[1,1]=cos(α)^2
        t[1,2]=cos(α)^2
        t[1,3]=-sin(α)*cos(α)
        t[2,1]=sin(α)^2
        t[2,2]=t[1,1]
        t[2,3]=-t[1,3]
        t[3,1]=2*sin(α)*cos(α)
        t[3,2]=-t[3,1]
        t[3,3]=cos(α)^2-sin(α)^2

        # Calcuate QT=C*T
        for i = 1:3
            for j = 1:3
                qt[i,j]=0.
                for k = 1:3
                    qt[i,j]=qt[i,j]+c[i,k]*t[k,j]
                end    
            end
        end        
        #println(qt)   
        # Calcuate \bar{q} = t*qt for [1,1] to [3,3]
        for i = 1:3
            for j = 1:3
                qbar[i,j]=0.
                for k = 1:3
                    qbar[i,j]=qbar[i,j]+t[k,i]*qt[k,j]
                end
            end
        end
        qbar[4,4]=g13*cos(α)^2+g23*sin(α)^2
        qbar[4,5]=(g13-g23)*cos(α)*sin(α)
        qbar[5,5]=g13*sin(α)^2+g23*cos(α)^2                 
        qbar[5,4]=qbar[4,5]

        # Compute a, b, d & ss matrix
        #println(qbar)
        for i = 1:3
            for j = 1:3
                a[i,j]=a[i,j]+qbar[i,j]*(zpos[l]-zpos[l+1])
            end
        end        

        for i = 1:3
            for j = 1:3
                b[i,j]=b[i,j]+qbar[i,j]*(zpos[l]^2-zpos[l+1]^2)
            end
        end

        for i = 1:3
            for j = 1:3
                d[i,j]=d[i,j]+qbar[i,j]*(zpos[l]^3-zpos[l+1]^3)
            end
        end

        ss[1,1]=ss[1,1]+qbar[4,4]*(zpos[l]-zpos[l+1])
        ss[1,2]=ss[1,2]+qbar[4,5]*(zpos[l]-zpos[l+1])
        ss[2,1]=ss[2,1]+qbar[5,4]*(zpos[l]-zpos[l+1])
        ss[2,2]=ss[2,2]+qbar[5,5]*(zpos[l]-zpos[l+1])

    end    
    #println(ss)
    for i = 1:3
        for j = 1:3
            cc[i,j]=a[i,j]
            cc[i,j+3]=b[i,j]/2.0
            cc[i+3,j]=b[i,j]/2.0
            cc[i+3,j+3]=d[i,j]/3.0
        end
    end
    
    cc[7,7]=5.0/6.0*ss[1,1]
    cc[7,8]=5.0/6.0*ss[1,2]        
    cc[8,7]=5.0/6.0*ss[2,1]
    cc[8,8]=5.0/6.0*ss[2,2]

    println(cc)
end
