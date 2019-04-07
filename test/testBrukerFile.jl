@testset "BrukerFile" begin

if !isdir("brukerfileCart")
  HTTP.open("GET", "http://media.tuhh.de/ibi/mrireco/brukerfileCart.zip") do http
    open("brukerfileCart.zip", "w") do file
        write(file, http)
    end
  end
  run(`unzip brukerfileCart.zip`)
end






end
