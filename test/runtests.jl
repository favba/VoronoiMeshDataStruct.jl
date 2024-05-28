using NCDatasets
using VoronoiMeshDataStruct
using TensorsLiteGeometry
using Test

@testset "VariableLengthStaticVector" begin

    @test_throws DomainError VariableLengthStaticVector{5}((1,0,2,3,5,0,0))

    @test_throws DomainError VariableLengthStaticVector{3}(VariableLengthStaticVector((-1,2,3,5,0,0),4))

    @test typeof(VariableLengthStaticVector((1.,2.,3.,4.,5.,0.,0.,0.),5)) === VariableLengthStaticVector{8,Float64}

    a = VariableLengthStaticVector((1.,2.,3.,4.,5.,0.,0.,0.),5)

    @test length(a) === 5
    @test size(a) === (5,)

    @test VoronoiMeshDataStruct.max_length(a) == 8
    @test VoronoiMeshDataStruct.max_length(typeof(a)) == 8

    @test collect(a) == [1., 2., 3., 4., 5.]

    @test_throws BoundsError getindex(a,6) 

    @test a[3] === 3.0
end

function compare_weights_trisk(m::Matrix,v::Vector{<:VariableLengthStaticVector})
    r = true
    for j in eachindex(v)
        vals = v[j]
        for i in eachindex(vals)
            r &= isapprox(vals[i],m[i,j],atol=1e-7)
        end
    end
    return r
end

@testset "Fields Creation Planar Mesh" begin
    mesh = NCDataset("mesh.nc") do f; VoronoiMesh(f) ;end
    xp = mesh.attributes[:x_period]::Float64
    yp = mesh.attributes[:y_period]::Float64

    @test all(x->(x[1]≈x[2]),  zip(compute_area_triangles(mesh),mesh.vertices.area))
    @test all(x->(x[1]≈x[2]),  zip(compute_area_cells(mesh),mesh.cells.area))
    @test all(x->(mapreduce(isapprox,&,x[1],x[2])),  zip(compute_kite_areas(mesh),mesh.vertices.kiteAreas))
    @test all(x->(x[1]≈x[2]),  zip(compute_dcEdge(mesh),mesh.edges.dc))
    @test all(x->(x[1]≈x[2]),  zip(compute_dvEdge(mesh),mesh.edges.dv))
    @test all(x->(x[1]≈x[2]),  zip(compute_angleEdge(mesh),mesh.edges.angle))
    @test all(x->(x[1]≈x[2]),  zip(compute_edge_position(mesh),mesh.edges.position))
    @test compare_weights_trisk(compute_weightsOnEdge_trisk(mesh),mesh.weightsOnEdge)

    @test length(find_obtuse_triangles(mesh)) == 0
end