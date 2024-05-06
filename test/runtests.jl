using NCDatasets
using VoronoiMeshDataStruct
using TensorsLiteGeometry
using Test

@testset "VariableLengthIndices" begin

    @test_throws DomainError VoronoiMeshDataStruct.VariableLengthIndices((1,0,2,3,5,0,0))

    @test_throws DomainError VoronoiMeshDataStruct.VariableLengthIndices((-1,2,3,5,0,0))

    @test typeof(VoronoiMeshDataStruct.VariableLengthIndices((1,2,3,4,5,0,0,0))) == VoronoiMeshDataStruct.VariableLengthIndices{8,Int}

    a = VoronoiMeshDataStruct.VariableLengthIndices((1,2,3,4,5,0,0,0))

    @test length(a) == 5
    @test VoronoiMeshDataStruct.max_length(a) == 8
    @test VoronoiMeshDataStruct.max_length(typeof(a)) == 8

    @test collect(a) == [1, 2, 3, 4, 5]

    @test_throws BoundsError getindex(a,6) 

    @test a[3] == 3

    @test VoronoiMeshDataStruct.VariableLengthIndices{6}((1,2,3)) === VoronoiMeshDataStruct.VariableLengthIndices((1,2,3,0,0,0))

    @test VoronoiMeshDataStruct.VariableLengthIndices{5}(a) === VoronoiMeshDataStruct.VariableLengthIndices((1,2,3,4,5))

    @test VoronoiMeshDataStruct.VariableLengthIndices{10}(a) === VoronoiMeshDataStruct.VariableLengthIndices((1,2,3,4,5,0,0,0,0,0))

    @test_throws DomainError VoronoiMeshDataStruct.VariableLengthIndices{3}(a)
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
    #@test all(x->(closest(x[2],x[1],xp,yp)≈x[2]),  zip(compute_edge_position(mesh),mesh.edges.position))
    @test all(x->(x[1]≈x[2]),  zip(compute_edge_position(mesh),mesh.edges.position))
    @test all(x->isapprox(x[1],x[2],atol=1e-7),  zip(compute_weightsOnEdge_trisk(mesh),mesh.weightsOnEdge))
end