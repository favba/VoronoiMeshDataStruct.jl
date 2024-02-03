using VoronoiMeshDataStruct
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