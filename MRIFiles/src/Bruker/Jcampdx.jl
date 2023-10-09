import Base: read, getindex, get, haskey

export JcampdxFile, findfirst_

findfirst_(A, v) = something(findfirst(isequal(v), A), 0)
function prevind_(A,i) 
  if i == 0
    return -1
  else
    return prevind(A,i)
  end
end

function repeatSubstringPattern(subLine::SubString{String})
  if(subLine[1]=='@')
    idx = findfirst_(subLine,'*')
    repN = subLine[2:idx-1]
    val = subLine[idx+2:end-1]
    
   return [val for _ in 1:parse(Int64,repN)]
  else
    return subLine
  end
end

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
      if line[1] == '#' && line[2] == '#' && line[3] == '$' 
        # the last command was not successfully parsed
        finishedReading = true
      end

      if finishedReading
        s = strip(line)
        if s[1] == '#' && s[2] == '#' && s[3] == '$'
          currentEntry += 1
          if maxEntries > 0 && currentEntry == maxEntries
            return file
          end

          i = findfirst_(s, '=')
          key = strip(s[4:prevind_(s,i)])

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
              if val[end] == ')' 
                file.dict[key] = val
              else
                remainingString = line[i+1:end]
                finishedReading = false
                currentKey = key
              end
            else
              j = findfirst_(val, ')')
              currentSizes = [parse(Int64,s) for s in split(val[2:prevind_(val,j)],",")]
              file.dict[key] = nothing
              currentKey = key
              currentIdx = 1
              finishedReading = false
            end
          end
        end
      else
         if line[1] == '<' && isnothing(remainingString)
           j = findfirst_(line, '>')
           file.dict[currentKey] = line[2:prevind_(line,j)]
           finishedReading = true
           tupleReading = false
         elseif line[1] == '(' || tupleReading  # skip these for the moment ...
           #finishedReading = true
           tupleReading = true

           if file.dict[currentKey] == nothing
             @debug "Will now allocate memory of size:" currentSizes
             file.dict[currentKey] = Array{Any}(undef, reverse(currentSizes)...)
           end

           totalLine = remainingString == nothing ? line[2:end] : remainingString*line

           parts = split(strip(totalLine), "(")

           for part in parts

             j = findfirst_(part, ')')
             if j != 0
               # Tuple read
               try
                 valsStr = split(strip(part[1:prevind_(part,j)]), ",")
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
         elseif remainingString !== nothing
          if line[end]==')'
            file.dict[currentKey] = remainingString*line[1:end]
            remainingString = nothing
            finishedReading = true
            currentKey = nothing
          else
            remainingString *= line
          end



         else
           #print(split(strip(line)," "))
           valsStr = split(strip(line), " ")
           vals=vcat(repeatSubstringPattern.(valsStr)...)
           
           if file.dict[currentKey] == nothing
             @debug "Will now allocate memory of size:" currentSizes
             file.dict[currentKey] = Array{eltype(vals)}(undef, reverse(currentSizes)...)
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
