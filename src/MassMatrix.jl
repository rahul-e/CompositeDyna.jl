function mass(xl::Vector{Float64}, yl::Vector{Float64}, th::Float64, rho::Float64)
    # xl and yl: x and y-coordinates of nodes in an element
    GP = Vector{Float64}(undef, 2)
    WG = Vector{Float64}(undef, 2)
    # Gaussian points
    GP=[-0.5773502691896, 0.5773502691896]
    WG=[1.0, 1.0]
    
    # X coordinates of the eight nodes in natural local coordinate system 
    # (refer Pg. 77, Fig. 5.14 in Finite Element Analysis by Bhavikatti)
    r = Array{Float64, 1}(undef, 8)
    r = [-1.0, 1.0, 1.0, -1.0, 0.0, 1.0, 0.0, -1.0]

    # Y coordinates of the eight nodes in natural local coordinate system 
    s = Array{Float64, 1}(undef, 8)
    s = [-1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 1.0, 0.0]

    # Initialize Jacobian matrix and its inverse 
    xjac=Array{Float64, 2}(undef, 2, 2)
    xjaci=Array{Float64, 2}(undef, 2, 2)
    for i=1:2
        for j=1:2
            xjac[i,j]=0.0
            xjaci[i,j]=0.0
        end
    end
    
    # Define vector for shape function
    sf = Array{Float64, 1}(undef, 8)
    sf = [i*0.0 for i in 1:8]
    # Define matrix for shape function derivatives
    # evaluated at Gaussian points 
    sfd = Array{Float64, 2}(undef, 2, 8)
    
    # Jacobian of shape function derivatives
    sfdj = Array{Float64, 2}(undef, 2, 8)
    
    # Initialize shape function and shape function derivatives
    for i=1:2
        for j=1:8
            sfd[i,j]=0.0
            sfdj[i,j]=0.0
        end
    end

    # Initialize determinant of Jacobian    
    Δjac::Float64=0.0
    da::Float64=0.0
    
    # Initialize element mass matrix
    em = Array{Float64, 2}(undef, 40, 40)
    for i=1:40
        for j=1:40
            em[i,j]=0.0
        end
    end
    
    # Define g
    g = Array{Float64, 2}(undef, 5, 40)
    sg = Array{Float64, 2}(undef, 5, 40)
    # Initialize miscellenous variables
    k1::Int64=0; k2::Int64=0; k3::Int64=0; k4::Int64=0; k5::Int64=0 

    #Define ss matrix
    ss = Array{Float64, 2}(undef, 5, 5)
    # Initialize ss matrix which contains terms that are functions of density and panel thickness   
    for i=1:5
        for j=1:5
            ss[i,j]=0.0
        end
    end 
    
    # Explicitly define diagonal terms (trace) of ss matrix  
    ss[1,1]=rho*th
    ss[2,2]=rho*th
    ss[3,3]=rho*th
    ss[4,4]=rho*th*th*th/12.0
    ss[5,5]=rho*th*th*th/12.0

    # Loop over to numerically integrate the terms in (SF[] × SS[]) using the four Gaussian-points whose weights are 1.0 
    for ix=1:2
        for iy=1:2
            # Calculate shape functions and their derivatives (refer Pg. 78 FEA by Bhavikatti)
            for i=1:8
                aa=(1.0+r[i]*GP[ix])
                bb=(1.0+s[i]*GP[iy])
                if(i<=4)
                    sf[i]=0.25*aa*bb*(aa+bb-3.0)
                    sfd[1,i]=0.25*bb*(2.0*aa+bb-3.0)*r[i]
                    sfd[2,i]=0.25*aa*(2.0*bb+aa-3.0)*s[i]
                else
                    aa=aa+(r[i]^2-1.0)*GP[ix]^2
                    bb=bb+(s[i]^2-1.0)*GP[iy]^2
                    sf[i]=0.5*aa*bb
                    sfd[1,i]=0.5*(r[i]+2.0*(r[i]^2-1.0)*GP[ix])*bb
                    sfd[2,i]=0.5*(s[i]+2.0*(s[i]^2-1.0)*GP[iy])*aa
                end
            end
            
                
            # Define elements in Jacobian matrix
            # Refer Pg. 227 Eqn. number 13.8 in FEA by Bhavikatti for 4-node quad element case.
            for i=1:8
                for j=1:2
                    xjac[j,1]=xjac[j,1]+xl[i]*sfd[j,i]
                    xjac[j,2]=xjac[j,2]+yl[i]*sfd[j,i]
                    #println(xl[i]," ",yl[i]," ",sfd[j,i])
                end
            end
            
            # Calculate determinant of Jacobian matrix
            Δjac=xjac[1,1]*xjac[2,2]-xjac[1,2]*xjac[2,1]
            
            # Define element in Inverse of Jacobian matrix
            xjaci[1,1]=xjac[2,2]/Δjac
            xjaci[1,2]=-1.0*xjac[1,2]/Δjac
            xjaci[2,1]=-1.0*xjac[2,1]/Δjac
            xjaci[2,2]=xjac[1,1]/Δjac

            da=Δjac*WG[ix]*WG[iy]
            #println(WG[ix], WG[ix], Δjac)
            for i=1:2
                for j=1:8
                    sfdj[i,j]=0.0
                end
            end

            for i=1:2
                for k=1:8
                    for j=1:2
                        sfdj[i,j]=sfdj[i,k]+xjaci[i,j]*sfd[j,k]
                    end
                end
            end


            for i=1:5
                for j=1:40
                    g[i,j]=0.0
                end
            end

            for i=1:8
                k1=5*(i-1)+1
                k2=k1+1
                k3=k2+1
                k4=k3+1
                k5=k4+1
                g[1,k1]=sf[i]
                g[2,k2]=sf[i]
                g[3,k3]=sf[i]
                g[4,k4]=sf[i]
                g[5,k5]=sf[i]
            end

            for i=1:5
                for j=1:40
                    sg[i,j]=0.0
                end
            end

            for i=1:5
                for j=1:40
                    for k=1:5
                        sg[i,j]=sg[i,j]+ss[i,k]*g[k,j]
                    end
                end
            end

            for i=1:40
                for j=1:40
                    for k=1:5
                        em[i,j]=em[i,j]+g[k,j]*sg[k,j]*da
                    end
                end
            end

        end
    end

    return em

end
