using Aqua, ExplicitImports

@testset "Aqua" begin
    Aqua.test_all(LambertWExp)
end

@testset "ExplicitImports" begin
    @test check_no_implicit_imports(LambertWExp) === nothing
    @test check_no_stale_explicit_imports(LambertWExp) === nothing
    @test check_all_explicit_imports_via_owners(LambertWExp) === nothing
    @test check_all_qualified_accesses_via_owners(LambertWExp) === nothing
    @test check_no_self_qualified_accesses(LambertWExp) === nothing
end
