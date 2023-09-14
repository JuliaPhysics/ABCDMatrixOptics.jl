using ABCDMatrixOptics, Documenter 

DocMeta.setdocmeta!(ABCDMatrixOptics, :DocTestSetup, :(using ABCDMatrixOptics); recursive=true)
makedocs(modules = [ABCDMatrixOptics], 
         sitename = "ABCDMatrixOptics.jl", 
         pages = Any[
            "ABCDMatrixOptics.jl" => "index.md",
            "Elements" => "elements.md",
            "Beams" => "beams.md",
            "Essential Functions" =>  "basics.md"
         ]
        )

deploydocs(repo = "github.com/JuliaPhysics/ABCDMatrixOptics.jl.git", devbranch="main")
