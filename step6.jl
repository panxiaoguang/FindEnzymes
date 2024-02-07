using BioSequences
using FASTX

function ParseExtraInfo(extrainfo::AbstractString)::Dict{String,Vector{NamedTuple{(:mismatch,:start,:seq),Tuple{Int64,Int64,LongSequence{DNAAlphabet{4}}}}}}
    result = Dict{String,Vector{NamedTuple{(:mismatch,:start,:seq),Tuple{Int64,Int64,LongSequence{DNAAlphabet{4}}}}}}()
    result["forward"] = Vector{NamedTuple{(:mismatch,:start,:seq),Tuple{Int64,Int64,LongSequence{DNAAlphabet{4}}}}}()
    result["reverse"] = Vector{NamedTuple{(:mismatch,:start,:seq),Tuple{Int64,Int64,LongSequence{DNAAlphabet{4}}}}}()
    if occursin(",", extrainfo)
        elements = split(extrainfo, ",")
        for ele in elements
            absloc,enzymeSeq,mis = split(ele,"|")
            if startswith(absloc,"+") # forward
                push!(result["forward"], (mismatch=parse(Int64, mis), start=parse(Int64, replace(absloc,"+"=>"")), seq=LongSequence{DNAAlphabet{4}}(enzymeSeq)))
            else # reverse
                push!(result["reverse"], (mismatch=parse(Int64, mis), start=parse(Int64, replace(absloc,"-"=>"")), seq=LongSequence{DNAAlphabet{4}}(enzymeSeq)))
            end
        end
    else
        absloc,enzymeSeq,mis = split(extrainfo,"|")
        if startswith(absloc,"+") # forward
            push!(result["forward"], (mismatch=parse(Int64, mis), start=parse(Int64, replace(absloc,"+"=>"")), seq=LongSequence{DNAAlphabet{4}}(enzymeSeq)))
        else # reverse
            push!(result["reverse"], (mismatch=parse(Int64, mis), start=parse(Int64, replace(absloc,"-"=>"")), seq=LongSequence{DNAAlphabet{4}}(enzymeSeq)))
        end
    end
    if result["forward"] == []
        result["forward"] = [(mismatch=0, start =0, seq=dna"NNN")]
    end
    if result["reverse"] == []
        result["reverse"] = [(mismatch=0, start =0, seq=dna"NNN")]
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
            name,seq = split(line, "\t")
            if name != site_name
                push!(enzymes, (name = name, seq = LongSequence{DNAAlphabet{4}}(seq)))
            end
        end
    end
    enzymes
end

function SearchRebase(dna::String, site_name::String, outF::String)
    dnas = open(FASTA.Reader, dna, index= string(dna, ".fai"))
    w = open(outF, "w")
    oldEnzymes = buildAllEnzymes("REBASE.txt", site_name)
    tmp_file = joinpath("Results/step4", site_name*".tsv")
    for line in eachline(tmp_file)
        id, chr, abpos, ref, alt,gene, sickType, extraInfo = split(line, '\t')
        oldExtraInfo = extraInfo
        dna = LongSequence{DNAAlphabet{4}}(FASTA.sequence(dnas[id]))
        _,_,pos,_,_,_,_,_,_,_,_ = split(FASTA.description(dnas[id]),"|")
        tmpEnzymeDict = ParseExtraInfo(extraInfo)
        pos = parse(Int64, pos)
        if tmpEnzymeDict["forward"] != [(mismatch=0, start =0, seq=dna"NNN")]
            upstream = dna[(pos-120):(pos-50-1)]
            downseq = dna[(pos+15+1):(pos+40)]
            for element in tmpEnzymeDict["forward"]
                result2 = find_restriction_site(downseq, element.seq)
                pushfirst!(oldEnzymes, (name = site_name, seq = element.seq))
                for (name,seq) in oldEnzymes
                    result1 = find_restriction_site(upstream, seq)
                    if result1 != Dict("forwad"=>[0], "reverse"=>[0]) && result2 != Dict("forwad"=>[0], "reverse"=>[0])
                        extraInfo = extraInfo * ":" * name * ","
                        break
                    end
                end
                popfirst!(oldEnzymes)
            end
        end
        if tmpEnzymeDict["reverse"] != [(mismatch=0, start =0, seq=dna"NNN")]
            downstream = dna[(pos+50+1):(pos+120)]
            upseq = dna[(pos-40):(pos-15-1)]
            for element in tmpEnzymeDict["reverse"]
                result2 = find_restriction_site(upseq, element.seq)
                pushfirst!(oldEnzymes, (name = site_name, seq = element.seq))
                for (name,seq) in oldEnzymes
                    result1 = find_restriction_site(downstream, seq)
                    if result1 != Dict("forwad"=>[0], "reverse"=>[0]) && result2 != Dict("forwad"=>[0], "reverse"=>[0])
                        extraInfo = extraInfo * ":" * name * ","
                        break
                    end
                end
                popfirst!(oldEnzymes)
            end
        end
        ## remove last ","
        if extraInfo != oldExtraInfo
            extraInfo = extraInfo[1:end-1]
            println(w, id, "\t", chr, "\t", abpos, "\t", ref, "\t", alt, "\t", gene, "\t", sickType, "\t", extraInfo)
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
        name, _ = split(line, "\t")
        name = convert(String,name)
        outname = joinpath(outpath,name * ".tsv")
        SearchRebase(inF, name, outname)
    end
end

main("REBASE.txt","snvs_modified.fasta","Results/step6")