function pasolv(sk::Array{Float64}, p::Array{Float64}, nds::Array{Int64}, nn::Int64, 
    neq1::Array{Int64}, nsky::Int64, inde::Int64)
    # PASOLV function from FEA - THeory and Programming, C. S. Krishnamoorthy Pg. 191
    # Direct solution routine using the algorithm based on Gauss elimination 
    # for static anaysis [K]{r}={P}
    n::Int64
    nn::Int64
    kn::Int64
    kl::Int64
    ku::Int64
    kh::Int64
    k::Int64
    ic::Int64
    klt::Int64

    for n=1:nn
        kn=nds[n]
        kl=kn+1
        ku=nds[n+1]-1
        kh=ku-kl
        if(kh>0)
            k=n-kh
            ic=0
            klt=ku
            for j=1:kh
                ic=ic+1
                klt=klt-1
                ki=nds[k]
                nd=nds[k+1]-ki-1
                if(nd>0) 
                    kk=min(ic,nd)
                    c=0.0
                    for l=1:kk
                        c=c+sk[ki+l]*sk[klt+l]
                    end
                    sk[klt]=sk[klt]-c    
                end
                k=k+1
            end
        elseif(kh==0)
            k=n
            b=0.0
            for kk=kl:ku
                k=k-1
                ki=nds[k]
                c=sk[kk]/sk[ki]
                b=b+c*sk[kk]
                sk[kk]=c
            end
            sk[kn]=sk[kn]-b
        else
            if(sk[kn]<=-0.001)
                println("n is %d", n)
                println("sk[%d] is %0.5f", kn, sk[kn])
                inde=1
                return
            end
        end
    end

    # Reduce right hand side load vector
    for n=1:nn
        kl=nds[n]+1
        ku=nds[n+1]-1

        if((ku-kl)>=0)
            k=n
            c=0.0
            for kk=kl:ku
                k=k-1
                c=c+sk[kk]*p[k]
            end
            p[n]=p[n]-c
        end
    end

    # Back substitution
    for n=1:nn
        k=nds[n]
        p[n]=p[n]/sk[k]
    end

    if(nn==1) 
        return
    end
    n=nn
    for l=2:nn
        kl=nds[n]+1
        ku=nds[n+1]-1
        if((ku-kl)>=0)
            k=n
        end
        for kk=kl:ku
            k=k-1
            p[k]=p[k]-sk[kk]*p[n]
        end
        n=n-1
    end
    return
end






