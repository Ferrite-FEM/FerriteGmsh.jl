using FerriteGmsh
using Ferrite

struct Elasticity
    λ::Float64
    μ::Float64
end

function constitutive_driver(ε::SymmetricTensor{2,dim}, material::Elasticity) where dim
    λ = material.λ
    μ = material.μ
    𝐈 = one(SymmetricTensor{2,dim})
    𝕀 = one(SymmetricTensor{4,dim})
    
    # Bestimmen sie hier σ und ℂ
    σ = λ * tr(ε) * 𝐈 + 2μ * ε  
    ℂ = λ * 𝐈 ⊗ 𝐈 + 2μ * 𝕀
    
    return σ, ℂ
end

function doassemble(cellvalues::CellValues, facevalues::FaceValues, dh::DofHandler,material::Elasticity) 
    traction = Dict("top" => Vec{2}((0,5)), "left" => Vec{2}((-5,0)))
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)

    K = create_sparsity_pattern(dh)
    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
    
    b = Vec{2}((0,0))
    
    @inbounds for (cellcount,cell) in enumerate(CellIterator(dh)) #für Element in Elemente
        fill!(Ke, 0)
        fill!(fe, 0)
        reinit!(cellvalues, cell)

        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)

            for i in 1:n_basefuncs
                δu  = shape_value(cellvalues, q_point, i)
                δε = shape_symmetric_gradient(cellvalues, q_point, i)
                fe[i] += (δu ⋅ b) * dΩ
                for j in 1:n_basefuncs
                    ε = shape_symmetric_gradient(cellvalues,q_point,j)
                    σ,ℂ = constitutive_driver(ε,material)
                    Ke[i, j] += (δɛ ⊡ σ) * dΩ
                end
            end
        end
        for face in 1:nfaces(cell)
            if ((cellcount, face) ∈ getfaceset(dh.grid, "left")) ||((cellcount, face) ∈ getfaceset(dh.grid, "top")) 
                t = Vec{2}((0.,0.))
                if ((cellcount, face) ∈ getfaceset(dh.grid, "left"))
                    t = traction["left"]
                elseif((cellcount, face) ∈ getfaceset(dh.grid, "top"))
                    t = traction["top"]
                end
                reinit!(facevalues, cell, face)
                for q_point in 1:getnquadpoints(facevalues)
                    dΓ = getdetJdV(facevalues, q_point)
                    for i in 1:n_basefuncs
                        δu = shape_value(facevalues, q_point, i)
                        fe[i] += (δu ⋅ t) * dΓ
                    end
                end
            end
        end
        assemble!(assembler, celldofs(cell), fe, Ke)
    end
    return K,f
end

function solve()
    λ = 42
    μ = 62.5
    material = Elasticity(λ,μ)
    ansatzfunktionen = Lagrange{2, RefTetrahedron, 1}()
    integration = QuadratureRule{2, RefTetrahedron}(2)
    elementvalues = CellVectorValues(integration, ansatzfunktionen)
    integration_rand = QuadratureRule{1, RefTetrahedron}(2)
    randvalues = FaceVectorValues(integration_rand, ansatzfunktionen)
    grid = saved_file_to_grid("test/quarter_plate.msh")
    addfaceset!(grid,"top", x->norm(x[2]) ≈ 50.0)
    addfaceset!(grid,"bottom", x->norm(x[2]) ≈ 0.0)
    addfaceset!(grid,"left", x->norm(x[1]) ≈ 0.0)
    addfaceset!(grid,"right", x->norm(x[1]) ≈ 50.0)
    
    dh = DofHandler(grid)
    #wir pushen die gesuchte Lösung mit Namen u und Dimension 2
    push!(dh, :u, 2)
    close!(dh)
    
    ch = ConstraintHandler(dh)
    ∂Ω = getfaceset(grid, "right")
    dbc = Dirichlet(:u, ∂Ω, (x, t) -> [0],[1])
    add!(ch, dbc)
    ∂Ω = getfaceset(grid, "bottom")
    dbc = Dirichlet(:u, ∂Ω, (x, t) -> [0],[2])
    add!(ch, dbc)
    close!(ch)
    update!(ch, 0.0)
    
    K, f = doassemble(elementvalues, randvalues, dh, material)
    apply!(K, f, ch)
    u = K \ f
    return u, dh
end

u, dh = solve()

vtk_grid("elasticity", dh) do vtk
    vtk_point_data(vtk, dh, u)
end
