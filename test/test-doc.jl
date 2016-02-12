module TestDoctest

## doctest
using MultipleTesting
using Base.Test

@unix ? begin
    using Lexicon
    dt = doctest(MultipleTesting)
    @test length(failed(dt)) == 0
end : println("! Skipping doctest on windows...")

end
