function geo(ac::Float64, bc::Float64, nx::Int64, ny::Int64, 
     kl::Int64, kr::Int64, kt::Int64, kb::Int64, 
     klp::Int64, krp::Int64, ktp::Int64, kbp::Int64, 
     #maxnode::Int64, maxelement::Int64, nelements::Int64,
     #nbig::Int64, nbigp::Int64, nnd::Array{Int64, 2}, nndp::Array{Int64, 2}, xco::Array{Float64, 2}, yco::Array{Float64, 2}, 
     #ndd::Array{Int64, 2}, nddp::Array{Int64, 2},
     #xnode::Int64, ynode::Int64, nhigh::Int64; nf=40, nfp=16
     )
          ndof=40
          nf=40 # 8 elements × 5 dofs (u, v, w, θ_x, θ_y )
          nfp=16 # 8 elements × 2 dofs (u, v)
          maxnode=500
          maxelement=50
           nelements=nx*ny
     
          nodemat = fill(0, (2*ny+1, 2*nx+1))
          xnode = fill(0.0, maxnode)
          ynode = fill(0.0, maxnode)
          #ndd = fill(0, (maxnode, 5))
          #nddp = fill(0, (maxnode, 2))
          nelmat = fill(0, (maxelement, 8))
          nnd = fill(0, (maxelement, nf))
          nndp = fill(0, (maxelement, nfp))
          xco = fill(0.0, (4, 8))
          yco = fill(0.0, (4, 8))
     
     
          count::Int64 = 0
          nhigh::Int64 = 0
          nbignode::Int64 = 0
          nbig::Int64 = 0
          nbigp::Int64 = 0
     
          # Element length along x is dx and y is dy  
          dx = ac/(1.0*nx)
          dy = bc/(1.0*ny)
     
          # Assign a node number to each physical node and store this in "nodemat"
          count = 1
          for i = 1:2:ny*2+1
               for j = 1:nx*2+1
                    nodemat[i,j]=count
                    count=count+1
               end
     
               if(i!=ny*2+1)
                    for j = 1:2:nx*2+1
                         nodemat[i+1,j]=count
                         count=count+1
                    end
                    
               end
          end
     
          #println(nodemat)
          
          nhigh=count-1
          for i = 1:ny*2+1
               for j = 1:nx*2+1
                    if(nodemat[i,j]!=0)
                         xnode[nodemat[i,j]]=(j-1)*dx/2.0
                         ynode[nodemat[i,j]]=(i-1)*dy/2.0
                             #println(nodemat[i,j])
                    end
               end
          end
          nbignode=nodemat[ny*2+1,nx*2+1]
     
          ndd = fill(1,(nbignode,5))
           nddp = fill(1,(nbignode,2))
          # Set boundary conditions (set value equal to zero corresponding to constrained DOF for every boundary node)
     
          # Left side boundary
          for i=1:ny*2+1
               if(klp==1)
                    nddp[nodemat[i,1],1]=0
               end
               
               if(klp==2)
                    nddp[nodemat[i,1],1]=0
                    nddp[nodemat[i,1],2]=0
               end
     
               if(kl==1)
                    ndd[nodemat[i,1],1]=0
                    ndd[nodemat[i,1],3]=0
                    ndd[nodemat[i,1],5]=0
               end
               
               if(kl==2)
                    ndd[nodemat[i,1],1]=0
                    ndd[nodemat[i,1],2]=0
                    ndd[nodemat[i,1],3]=0
                    ndd[nodemat[i,1],4]=0
                    ndd[nodemat[i,1],5]=0
               end     
               
          end          
     
          # Right side boundary
          for i=1:ny*2+1
               if(krp==1)
                    nddp[nodemat[i,nx*2+1],1]=0
               end
               
               if(krp==2)
                    nddp[nodemat[i,nx*2+1],1]=0
                    nddp[nodemat[i,nx*2+1],2]=0
               end
     
               if(kr==1)
                    ndd[nodemat[i,nx*2+1],1]=0
                    ndd[nodemat[i,nx*2+1],3]=0
                    ndd[nodemat[i,nx*2+1],5]=0
               end
               
               if(kr==2)
                    ndd[nodemat[i,nx*2+1],1]=0
                    ndd[nodemat[i,nx*2+1],2]=0
                    ndd[nodemat[i,nx*2+1],3]=0
                    ndd[nodemat[i,nx*2+1],4]=0
                    ndd[nodemat[i,nx*2+1],5]=0
               end     
               
          end
          
          # Bottom side boundary
          for i=1:nx*2+1
               if(kbp==1)
                    nddp[nodemat[1,i],1]=0
               end
               
               if(kbp==2)
                    nddp[nodemat[1,i],1]=0
                    nddp[nodemat[1,i],2]=0
               end
     
               if(kb==1)
                    ndd[nodemat[1,i],2]=0
                    ndd[nodemat[1,i],3]=0
                    ndd[nodemat[1,i],4]=0
               end
               
               if(kb==2)
                    ndd[nodemat[1,i],1]=0
                    ndd[nodemat[1,i],2]=0
                    ndd[nodemat[1,i],3]=0
                    ndd[nodemat[1,i],4]=0
                    ndd[nodemat[1,i],5]=0
               end     
               
          end
     
          # Top side boundary
          for i=1:nx*2+1
               if(ktp==1)
                    nddp[nodemat[ny*2+1,i],1]=0
               end
               
               if(ktp==2)
                    nddp[nodemat[ny*2+1,i],1]=0
                    nddp[nodemat[ny*2+1,i],2]=0
               end
     
               if(kt==1)
                    ndd[nodemat[ny*2+1,i],2]=0
                    ndd[nodemat[ny*2+1,i],3]=0
                    ndd[nodemat[ny*2+1,i],4]=0
               end
               
               if(kt==2)
                    ndd[nodemat[ny*2+1,i],1]=0
                    ndd[nodemat[ny*2+1,i],2]=0
                    ndd[nodemat[ny*2+1,i],3]=0
                    ndd[nodemat[ny*2+1,i],4]=0
                    ndd[nodemat[ny*2+1,i],5]=0
               end     
               
          end
     
          # Replace placeholder (value one) in ndd and nddp by equation numbers
     
          nbig=0
          nbigp=0
          
          for i=1:nbignode
               for j=1:5
                    if(ndd[i,j]!=0)
                         nbig=nbig+1
                         ndd[i,j]=nbig
                    end
               end
               
               for j=1:2
                    if(nddp[i,j]!=0)
                         nbigp=nbigp+1
                         nddp[i,j]=nbigp
                    end
               end      
          end
     
          # Create "nelmat" which is a map for each element's local node number to the global 
          # node number stored in "nodemat". Nodes of an element are numbered in the ordered
          # [1, 2, 3, 4] for corner nodes beginning from bottom-left counting anti-clockwise
          # and [5, 6, 7, 8] for midside nodes starting from bottom side of element.  
     
          for i=1:ny
               for j=1:nx
                    num=(i-1)*nx+j
                    nelmat[num,1]=nodemat[(i-1)*2+1, (j-1)*2+1]
                    nelmat[num,2]=nodemat[(i-1)*2+1, (j-1)*2+3]
                    nelmat[num,3]=nodemat[(i-1)*2+3, (j-1)*2+3]
                    nelmat[num,4]=nodemat[(i-1)*2+3, (j-1)*2+1]
                    nelmat[num,5]=nodemat[(i-1)*2+1, (j-1)*2+2]
                    nelmat[num,6]=nodemat[(i-1)*2+2, (j-1)*2+3]               
                    nelmat[num,7]=nodemat[(i-1)*2+3, (j-1)*2+2]
                    nelmat[num,8]=nodemat[(i-1)*2+2, (j-1)*2+1]
               end
          end          
          
          # Total number of elements
          nelements=nx*ny
          for i=1:nelements
               for j=1:8
                    for k=1:5
                         num=(j-1)*5+k
                         nnd[i,num]=ndd[nelmat[i,j],k]
                    end
     
                    for k=1:2
                         num=(j-1)*2+k
                         nndp[i,num]=nddp[nelmat[i,j],k]
                    end
               end
          end
     
           # Local (x,y) coordinates 
     
          for i=1:4
               
               xco[i,1]=0.0
               xco[i,2]=dx
               xco[i,3]=dx
               xco[i,4]=0.0
               xco[i,5]=dx/2.0
               xco[i,6]=dx
               xco[i,7]=dx/2.0
               xco[i,8]=0.0
     
               yco[i,1]=0.0
               yco[i,2]=0.0
               yco[i,3]=dy
               yco[i,4]=dy
               yco[i,5]=0.0
               yco[i,6]=dy/2.0
               yco[i,7]=dy
               yco[i,8]=dy/2.0
     
          end
          
           xl = fill(0.0, 8)
           yl = fill(0.0, 8)
          
          for i=1:8
               xl[i]=xco[1,i]
               yl[i]=yco[1,i]
          end
          
          # Call BendingStiff.jl
          # Call MassMatrix.jl
     
          # Initiate skyline matrix
          nsky::Int64=0
          
           cht = fill(0, nbig)
           nds = fill(0, nbig+1)
          
          #mband::Int64=0
          lnum::Int64=0
          # CHT = Column Height
          nd = fill(0, ndof)
          
          for lnum=1:nelements
               for i=1:40
                    nd[i]=nnd[lnum,i]
                 end
               cht = ColumnH(cht, nd, 40, nbig)
          end
          for lnum=1:nelements
               #nbig1=nbig+1
               nds, nsky = CadNum(cht, nbig, nsky)
          end
     
          # Initialize global stiffness and mass matrix and
           gk = fill(0.0, nsky)
           gm = fill(0.0, nsky)
          
          sqkg = fill(0.0, (nbig, nbig))
           sqk = fill(0.0, (nbig, nbig))
                
           # Call the funtion to generate element mass matrix
           em = mass(xl,yl,th)
           ek = bstif(cc, xl, yl, 1.0E10, 1.0E10)
           
           # Loop over elements to assemble
          for lnum=1:nelements
               for i=1:40
                    nd[i]=nnd[lnum,i]
                 end
                 
               for i=1:40
                    for j=1:40
                         if(nd[i]!=0 && nd[j]!=0)
                              sqk[nd[i],nd[j]]=sqk[nd[i],nd[j]]+em[i,j]
                         end
                    end
               end
     
               nxl=mod(lnum,nx)
               nyl=((lnum-1)/nx)+1
     
               if(nxl==0)
                    nxl=nx
               end
     
               #if((nxl<=idxe) && (nxl>=idxs) && (nyl<=idye) && (nyl>=idys))
               #     passem(GK,DEK,NDS,ND,40,NBIG1,NSKY,40)
               #else
               gk=passem(gk,ek,nds,nd,40,nbig+1,nsky,40)
               #end
     
               gm=passem(gm,em,nds,nd,40,nbig+1,nsky,40)
          end
     
          nn::Int64=nbig
          nwk::Int64=nsky
          println("Number of elements in the stiffness matrix is ", nsky)
           println("Number of equations in the global system is ", nbig)
          println("Computing free vibration characteristics")
     
          return gk, gm, nbig, nsky, cht
     end  