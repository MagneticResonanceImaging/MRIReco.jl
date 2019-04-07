import Base: read, getindex, get, haskey

export JcampdxFile, findfirst_

findfirst_(A, v) = something(findfirst(isequal(v), A), 0)

const JCVAL = Union{AbstractString,Number,Bool,Array,Tuple,Nothing}
const HTSS = Dict{AbstractString,JCVAL}

mutable struct JcampdxFile
  dict::HTSS
end

JcampdxFile() = JcampdxFile(HTSS())

getindex(file::JcampdxFile, key::AbstractString) = file.dict[key]
haskey(file::JcampdxFile, key::AbstractString) = haskey(file.dict, key)
get(file::JcampdxFile, key::AbstractString, default) = get(file.dict, key, default)

function read(file::JcampdxFile, stream::IO, keylist::Vector=String[]; maxEntries=-1)
    finishedReading = true
    currentKey = nothing
    currentIdx = 1
    currentSizes = nothing
    currentEntry = 0
    remainingString = nothing
    tupleReading = false

    skipKeys = ["VisuAcqFrameTime"]

    for line in eachline(stream)
      if finishedReading
        s = strip(line)
        if s[1] == '#' && s[2] == '#' && s[3] == '$'
          currentEntry += 1
          if maxEntries > 0 && currentEntry == maxEntries
            return file
          end

          i = findfirst_(s, '=')
          key = strip(s[4:i-1])

          # Small HACK
          if in(key, skipKeys)
            return file
          end

          if !isempty(keylist) && findfirst_(keylist,key) == 0
            continue
          end

          val = strip(s[i+1:end])
          if val[1] != '('
            file.dict[key] = val
          else
            if val[2] != ' '
              file.dict[key] = val
            else
              j = findfirst_(val, ')')
              currentSizes = [parse(Int64,s) for s in split(val[2:j-1],",")]
              file.dict[key] = nothing
              currentKey = key
              currentIdx = 1
              finishedReading = false
            end
          end
        end
      else
         if line[1] == '<'
           j = findfirst_(line, '>')
           file.dict[currentKey] = line[2:j-1]
           finishedReading = true
           tupleReading = false
         elseif line[1] == '(' || tupleReading  # skip these for the moment ...
           #finishedReading = true
           tupleReading = true

           if file.dict[currentKey] == nothing
             @debug "Will now allocate memory of size:" currentSizes
             file.dict[currentKey] = Array{Any}(undef, currentSizes...)
           end

           totalLine = remainingString == nothing ? line[2:end] : remainingString*line

           parts = split(strip(totalLine), "(")

           for part in parts

             j = findfirst_(part, ')')
             if j != 0
               # Tuple read
               try
                 valsStr = split(strip(part[1:j-1]), ",")
                 vals = Any[]
                 for valStr in valsStr
                   #try
                   #  push!(vals, int(valStr))
                   #catch
                     try
                       push!(vals, parse(Float64,valStr) )
                     catch
                       push!(vals, valStr)
                     end
                   #end
                 end

                 file.dict[currentKey][currentIdx] = vals
                 currentIdx += 1
                 if currentIdx-1 >= length(file.dict[currentKey])
                   finishedReading = true
                   tupleReading = false
                 end
               catch
                 @debug "" currentIdx vals currentSizes size(file.dict[currentKey]) length(keys(file.dict)) currentKey
                 #rethrow()
                 finishedReading = true
                 tupleReading = false
               end

               remainingString = nothing
             else
               remainingString = part
             end
           end




         else
           #print(split(strip(line)," "))
           valsStr = split(strip(line), " ")
           vals = nothing
           try
             vals = parse(Int64,valsStr)
           catch
             try
               vals = parse(Float64,valsStr)
             catch
               vals = valsStr
             end
           end

           if file.dict[currentKey] == nothing
             @debug "Will now allocate memory of size:" currentSizes
             file.dict[currentKey] = Array{eltype(vals)}(undef, currentSizes...)
           end

           try
             file.dict[currentKey][currentIdx:currentIdx+length(vals)-1] = vals
             currentIdx += length(vals)
             if currentIdx-1 >= length(file.dict[currentKey])
               finishedReading = true
               tupleReading = false
             end
           catch
             @debug "" currentIdx vals currentSizes size(file.dict[currentKey])
             rethrow()
           end
         end
      end
    end
   file
end

function read(file::JcampdxFile, filename::AbstractString,
              keylist::Vector=String[]; maxEntries=-1)
    open(filename) do f
        read(file, f, keylist, maxEntries=maxEntries)
    end
    file
end
