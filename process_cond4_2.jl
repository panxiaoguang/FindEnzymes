using BioSequences
using FASTX

function ParseExtraInfo(extrainfo::AbstractString)::Dict{String,Vector{NamedTuple{(:orient,:start,:seq),Tuple{Int64,Int64,LongSequence{DNAAlphabet{4}}}}}}
    result = Dict{String,Vector{NamedTuple{(:orient,:start,:seq),Tuple{Int64,Int64,LongSequence{DNAAlphabet{4}}}}}}()
    result["forward"] = Vector{NamedTuple{(:orient,:start,:seq),Tuple{Int64,Int64,LongSequence{DNAAlphabet{4}}}}}()
    result["reverse"] = Vector{NamedTuple{(:orient,:start,:seq),Tuple{Int64,Int64,LongSequence{DNAAlphabet{4}}}}}()
    if occursin(",", extrainfo)
        elements = split(extrainfo, ",")
        for ele in elements
            absloc,enzymeSeq,orient = split(ele,"|")
            if startswith(absloc,"+") # forward
                push!(result["forward"], (orient=parse(Int64, orient), start=parse(Int64, replace(absloc,"+"=>"")), seq=LongSequence{DNAAlphabet{4}}(enzymeSeq)))
            else # reverse
                push!(result["reverse"], (orient=parse(Int64, orient), start=parse(Int64, replace(absloc,"-"=>"")), seq=LongSequence{DNAAlphabet{4}}(enzymeSeq)))
            end
        end
    else
        absloc,enzymeSeq,orient = split(extrainfo,"|")
        if startswith(absloc,"+") # forward
            push!(result["forward"], (orient=parse(Int64, orient), start=parse(Int64, replace(absloc,"+"=>"")), seq=LongSequence{DNAAlphabet{4}}(enzymeSeq)))
        else # reverse
            push!(result["reverse"], (orient=parse(Int64, orient), start=parse(Int64, replace(absloc,"-"=>"")), seq=LongSequence{DNAAlphabet{4}}(enzymeSeq)))
        end
    end
    if result["forward"] == []
        result["forward"] = [(orient=0, start =0, seq=dna"NNN")]
    end
    if result["reverse"] == []
        result["reverse"] = [(orient=0, start =0, seq=dna"NNN")]
    end
    result
end

function find_restriction_site(dna::LongSequence{DNAAlphabet{4}}, site::LongSequence{DNAAlphabet{4}})::Dict{String,Vector{Int64}}
    query = ExactSearchQuery(site, iscompatible)
    reverse_dna = reverse_complement(dna)
    forward_check = occursin(query, dna)
    reverse_check = occursin(query, reverse_dna)
    CheckResult = Dict("forwad"=>Vector{Int64}(), "reverse"=>Vector{Int64}())
    if forward_check && !reverse_check
        locs = findall(query, dna)
        startSite = [loc.start for loc in locs]
        CheckResult["forwad"] = startSite
        CheckResult["reverse"] = [0]
    elseif !forward_check && reverse_check
        locs_reverse = findall(query, reverse_dna)
        startSite_reverse = [loc.start for loc in locs_reverse]
        CheckResult["forwad"] = [0]
        CheckResult["reverse"] = startSite_reverse
    elseif forward_check && reverse_check
        locs = findall(query, dna)
        startSite = [loc.start for loc in locs]
        locs_reverse = findall(query, reverse_dna)
        startSite_reverse = [loc.start for loc in locs_reverse]
        CheckResult["forwad"] = startSite
        CheckResult["reverse"] = startSite_reverse
    else
        CheckResult["forwad"] = [0]
        CheckResult["reverse"] = [0]
    end
    CheckResult
end


function buildAllEnzymes(enzymeFile::String,site_name::String)::Vector{NamedTuple{(:name,:seq),Tuple{String,LongSequence{DNAAlphabet{4}}}}}
    enzymes = Vector{NamedTuple{(:name,:seq),Tuple{String,LongSequence{DNAAlphabet{4}}}}}()
    open(enzymeFile,"r") do file
        for line in eachline(file)
            if startswith(line,"Name")
                continue
            end
            name,seq = split(line, "\t")
            name = convert(String,name)
            seq = LongSequence{DNAAlphabet{4}}(seq)
            if name != site_name
                push!(enzymes, (name = name, seq = seq))
            end
        end
    end
    enzymes
end


function SearchRebase(dna::String, site_name::String, outF::String)
    dnas = open(FASTA.Reader, dna, index= string(dna, ".fai"))
    w = open(outF, "w")
    oldEnzymes = buildAllEnzymes("REs.tsv", site_name)
    tmp_file = joinpath("firstStep/cond4_1", site_name*".tsv")
    for line in eachline(tmp_file)
        id, chr, abpos, ref, alt,gene, hgvsg, extraInfo = split(line, '\t')
        oldExtraInfo = extraInfo
        dna = LongSequence{DNAAlphabet{4}}(FASTA.sequence(dnas[id]))
        des = FASTA.description(dnas[id])
        pos = parse(Int64,split(des," ")[2])
        tmpEnzymeDict = ParseExtraInfo(extraInfo)
        upstream = dna[(pos-120):(pos-50-1)]
        downstream = dna[(pos+50+1):(pos+120)]
        if tmpEnzymeDict["forward"] != [(orient=0, start =0, seq=dna"NNN")]
            for element in tmpEnzymeDict["forward"]
                pushfirst!(oldEnzymes, (name = site_name, seq = element.seq))
                if element.orient == 1
                    for (name,seq) in oldEnzymes
                        result = find_restriction_site(upstream, seq)
                        if result != Dict("forwad"=>[0], "reverse"=>[0])
                            extraInfo = extraInfo * ":" * name * ","
                            break
                        end
                    end
                else
                    for (name,seq) in oldEnzymes
                        result = find_restriction_site(downstream, seq)
                        if result != Dict("forwad"=>[0], "reverse"=>[0])
                            extraInfo = extraInfo * ":" * name * ","
                            break
                        end
                    end
                end
                popfirst!(oldEnzymes)
            end
        end
        if tmpEnzymeDict["reverse"] != [(orient=0, start =0, seq=dna"NNN")]
            for element in tmpEnzymeDict["reverse"]
                pushfirst!(oldEnzymes, (name = site_name, seq = element.seq))
                if element.orient == 1
                    for (name,seq) in oldEnzymes
                        result = find_restriction_site(upstream, seq)
                        if result != Dict("forwad"=>[0], "reverse"=>[0])
                            extraInfo = extraInfo * ":" * name * ","
                            break
                        end
                    end
                else
                    for (name,seq) in oldEnzymes
                        result = find_restriction_site(downstream, seq)
                        if result != Dict("forwad"=>[0], "reverse"=>[0])
                            extraInfo = extraInfo * ":" * name * ","
                            break
                        end
                    end
                end
                popfirst!(oldEnzymes)
            end
        end
        ## remove last ","
        if extraInfo != oldExtraInfo
            extraInfo = extraInfo[1:end-1]
            println(w, id, "\t", chr, "\t", abpos, "\t", ref, "\t", alt, "\t", gene, "\t", hgvsg, "\t", extraInfo)
        end
    end
    close(dnas)
    close(w)
end

### Test

#SearchRebase("snvs_modified.fasta","AarI")

function main(rebase::String,inF::String,outpath::String)
    if !ispath(outpath)
        mkpath(outpath)
    end
    for line in eachline(rebase)
        if startswith(line,"Name")
            continue
        end
        name, _ = split(line, "\t")
        name = convert(String,name)
        outname = joinpath(outpath,name * ".tsv")
        SearchRebase(inF, name, outname)
    end
end

main("REs.tsv","Raw_sequences.fasta","firstStep/cond4_2")


