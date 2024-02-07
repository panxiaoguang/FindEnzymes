using BioSequences
using FASTX

function find_first_mismatch(longseq::LongSequence{DNAAlphabet{4}}, shortseq::LongSequence{DNAAlphabet{4}}, start::Int64, center::Int64)::Int64
    dist = 0
    if length(shortseq)+start-1 > length(longseq)
        dist = 0
    else
        truncted_seq = longseq[start:(length(shortseq)+start-1)]
        for i in eachindex(truncted_seq)
            if !iscompatible(truncted_seq[i],shortseq[i])
                dist = i - (center-start + 1) - 1
                break
            end
        end
    end
    dist
end
function backtrace_sequence(longseq::LongSequence{DNAAlphabet{4}}, shortseq::LongSequence{DNAAlphabet{4}}, start::Int64)::LongSequence{DNAAlphabet{4}}
    truncted_seq = longseq[start:(length(shortseq)+start-1)]
    newseq = copy(shortseq)
    for i in eachindex(truncted_seq)
        if !iscompatible(truncted_seq[i],shortseq[i])
            newseq[i] = truncted_seq[i]
        end
    end
    newseq
end

function match_sequence(seq::LongSequence{DNAAlphabet{4}},pattern::LongSequence{DNAAlphabet{4}},max_mismatch::Int64)
    query = ApproximateSearchQuery(pattern, iscompatible)
    center = convert(Int64,(length(seq)-1)/2+1)
    result = Dict{String,Vector{NamedTuple{(:mismatch,:start,:seq),Tuple{Int64,Int64,LongSequence{DNAAlphabet{4}}}}}}()
    result["forward"] = Vector{NamedTuple{(:mismatch,:start,:seq),Tuple{Int64,Int64,LongSequence{DNAAlphabet{4}}}}}()
    result["reverse"] = Vector{NamedTuple{(:mismatch,:start,:seq),Tuple{Int64,Int64,LongSequence{DNAAlphabet{4}}}}}()
    #result = Vector{NamedTuple{(:mismatch,:start,:seq),Tuple{Int64,Int64,LongSequence{DNAAlphabet{4}}}}}()
    revSeq = reverse_complement(seq)
    for mis in 1:max_mismatch
        start = 1
        while true
            pipei = findnext(query,mis,seq,start)
            if pipei == nothing
                break
            elseif pipei.stop-pipei.start+1 == length(pattern) && find_first_mismatch(seq,pattern,pipei.start,center) >= 3
                push!(result["forward"],(mismatch=mis, start =pipei.start, seq=backtrace_sequence(seq,pattern,pipei.start)))
            end
            start = pipei.start + 1
        end
    end
    if result["forward"] == []
        result["forward"] = [(mismatch=0, start =0, seq=dna"NNN")]
    end
    for mis in 1:max_mismatch
        start = 1
        while true
            pipei = findnext(query,mis,revSeq,start)
            if pipei == nothing
                break
            elseif pipei.stop-pipei.start+1 == length(pattern) && find_first_mismatch(revSeq,pattern,pipei.start,center) >= 3
                push!(result["reverse"],(mismatch=mis, start =pipei.start, seq=backtrace_sequence(revSeq,pattern,pipei.start)))
            end
            start = pipei.start + 1
        end
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

function get_step1_IDs(site_name::String)::Vector{String}
    first_step_file = joinpath("Results/step1", site_name*".tsv")
    ids = Vector{String}()
    for line in eachline(first_step_file)
        id, _, _, _, _,_, _, _ = split(line, '\t')
        id = convert(String, id)
        push!(ids, id)
    end
    ids
end


function printTabule(outf::IO,matched::Dict{String,Vector{NamedTuple{(:mismatch,:start,:seq),Tuple{Int64,Int64,LongSequence{DNAAlphabet{4}}}}}},mutate_dna::LongSequence{DNAAlphabet{4}},id::String,chr::String,abpos::Int64,ref::String,alt::String,gene::String,sickType::String,site_length::Int64)
    if matched != Dict("forward"=>[(mismatch=0, start =0, seq=dna"NNN")],"reverse"=>[(mismatch=0, start =0, seq=dna"NNN")])
        extraInfo = ""
        if matched["forward"] != [(mismatch=0, start =0, seq=dna"NNN")]
            for tp in matched["forward"]
                if find_restriction_site(mutate_dna,tp.seq) == Dict("forwad"=>[0], "reverse"=>[0])
                    absloc = abpos-site_length+(tp.start-1)
                    extraInfo = extraInfo * "+" * string(absloc) *"|"* string(tp.seq)*"|"*string(tp.mismatch)*","
                end
            end
        end
        if matched["reverse"] != [(mismatch=0, start =0, seq=dna"NNN")]
            for tp in matched["reverse"]
                if find_restriction_site(mutate_dna,tp.seq) == Dict("forwad"=>[0], "reverse"=>[0])
                    absloc = abpos+site_length-(tp.start-1)
                    extraInfo = extraInfo * "-" * string(absloc) *"|"* string(tp.seq)*"|"*string(tp.mismatch)*","
                end
            end
        end
        if extraInfo != ""
            extraInfo = extraInfo[1:end-1]
            println(outf,id,"\t",chr,"\t",abpos,"\t",ref,"\t",alt,"\t",gene,"\t",sickType,"\t",extraInfo)
        end
    end
end


function SearchRebase(dna::String, site_name::String, site::LongSequence{DNAAlphabet{4}},outF::String)
    dnas = open(FASTA.Reader, dna, index= string(dna, ".fai"))
    filterID = get_step1_IDs(site_name)
    w = open(outF, "w")
    ## read the file from step1
    for record in dnas
        if (FASTA.identifier(record) in filterID)
            continue
        else
            dna = LongSequence{DNAAlphabet{4}}(FASTA.sequence(record))
            id = FASTA.identifier(record)
            id = convert(String, id)
            name = FASTA.description(record)
            chr,abpos,pos,ref,alt,_,sick,_,_,gene,_ = split(name,"|")
            chr = replace(chr,r"^\d+\s"=>"")
            sick = convert(String, sick)
            chr = convert(String, chr)
            abpos = parse(Int64, abpos)
            ref = convert(String, ref)
            alt = convert(String, alt)
            gene = convert(String, gene)
            if alt=="."
                continue
            else
                pos = parse(Int, pos)
                site_length = length(site)
                dna = dna[(pos-site_length):(pos+site_length)]
                mutate_dna = copy(dna)
                center_coord = convert(Int64,(length(dna)-1) /2 +1)
                mutate_dna[center_coord] = convert(DNA, alt[1])
                matched = match_sequence(dna, site, 3)
                printTabule(w, matched,mutate_dna,id,chr,abpos,ref,alt,gene,sick,site_length)
            end
        end
    end
    close(dnas)
    close(w)
end

## test 

#SearchRebase("snvs_modified.fasta","FalI",dna"AAGNNNNNCTT","test4.tsv")

function main(rebase::String,inF::String,outpath::String)
    if !ispath(outpath)
        mkpath(outpath)
    end
    for line in eachline(rebase)
        name, site = split(line, "\t")
        name = convert(String,name)
        site = LongSequence{DNAAlphabet{4}}(site)
        outname = joinpath(outpath,name * ".tsv")
        SearchRebase(inF, name, site, outname)
    end
end

main("REBASE.txt","snvs_modified.fasta","Results/step4")