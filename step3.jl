using BioSequences
using FASTX


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

function isMatch(locs::Dict{String,Vector{Int64}})::Bool
    if locs != Dict("forwad"=>[0], "reverse"=>[0])
        return true
    else
        return false
    end
end

function printTabule(outf::IO, uplocs::Dict{String,Vector{Int64}},downlocs::Dict{String,Vector{Int64}},uplocs2::Dict{String,Vector{Int64}},downlocs2::Dict{String,Vector{Int64}},id::String,sickType::String,chr::String,abpos::Int64,ref::String,alt::String,gene::String,oldInfo::String)
    if (isMatch(uplocs) && isMatch(downlocs)) || (isMatch(uplocs2) && isMatch(downlocs2))
        extraInfo = ""
        if uplocs["forwad"] != [0]
            for loc in uplocs["forwad"]
                ## should calculate the absolute position
                absloc = abpos-120+(loc-1)
                extraInfo = extraInfo* "up10|+" * string(absloc) * ","
            end
        end
        if uplocs["reverse"] != [0]
            for loc in uplocs["reverse"]
                absloc = abpos-(10-1)-(loc-1)
                extraInfo = extraInfo* "up10|-" * string(absloc) * ","
            end
        end
        if downlocs["forwad"] != [0]
            for loc in downlocs["forwad"]
                ## should calculate the absolute position
                absloc = abpos+(50+1)+(loc-1)
                extraInfo = extraInfo* "down50|+" * string(absloc) * ","
            end
        end
        if downlocs["reverse"] != [0]
            for loc in downlocs["reverse"]
                absloc = abpos+120-(loc-1)
                extraInfo = extraInfo* "down50|-" * string(absloc) * ","
            end
        end
        if uplocs2["forwad"] != [0]
            for loc in uplocs2["forwad"]
                ## should calculate the absolute position
                absloc = abpos-120+(loc-1)
                extraInfo = extraInfo* "up50|+" * string(absloc) * ","
            end
        end
        if uplocs2["reverse"] != [0]
            for loc in uplocs2["reverse"]
                absloc = abpos-(50-1)-(loc-1)
                extraInfo = extraInfo* "up50|-" * string(absloc) * ","
            end
        end
        if downlocs2["forwad"] != [0]
            for loc in downlocs2["forwad"]
                ## should calculate the absolute position
                absloc = abpos+(10+1)+(loc-1)
                extraInfo = extraInfo* "down10|+" * string(absloc) * ","
            end
        end
        if downlocs2["reverse"] != [0]
            for loc in downlocs2["reverse"]
                absloc = abpos+120-(loc-1)
                extraInfo = extraInfo* "down10|-" * string(absloc) * ","
            end
        end
        extraInfo = extraInfo[1:end-1]
        println(outf, id, "\t", chr, "\t", abpos, "\t", ref, "\t", alt, "\t", gene, "\t", sickType, "\t", oldInfo, "\t", extraInfo)
    end
end

function SearchRebase(dna::String, site_name::String, site::LongSequence{DNAAlphabet{4}},outF::String)
    dnas = open(FASTA.Reader, dna, index= string(dna, ".fai"))
    w = open(outF, "w")
    ## read the file from step1
    first_step_file = joinpath("Results/step1", site_name*".tsv")
    for line in eachline(first_step_file)
        id, chr, abpos, ref, alt,gene, sickType, extraInfo = split(line, '\t')
        id = convert(String, id)
        chr = convert(String, chr)
        abpos = parse(Int64, abpos)
        ref = convert(String, ref)
        alt = convert(String, alt)
        gene = convert(String, gene)
        sickType = convert(String, sickType)
        extraInfo = convert(String, extraInfo)
        ## get the sequence from the id and dna
        dna = LongSequence{DNAAlphabet{4}}(FASTA.sequence(dnas[id]))
        _,_,pos,_,_,_,_,_,_,_,_ = split(FASTA.description(dnas[id]),"|")
        ## check from -120bp to 10bp
        pos = parse(Int64, pos)
        upstream_1 = dna[(pos-120):(pos-10-1)]
        up_locs_1 = find_restriction_site(upstream_1, site)
        ## check from 10bp to 120bp
        downstream_1 = dna[(pos+50+1):(pos+120)]
        down_locs_1 = find_restriction_site(downstream_1, site)
        upstream_2 = dna[(pos-120):(pos-50-1)]
        up_locs_2 = find_restriction_site(upstream_2, site)
        ## check from 10bp to 120bp
        downstream_2 = dna[(pos+10+1):(pos+120)]
        down_locs_2 = find_restriction_site(downstream_2, site)
        printTabule(w, up_locs_1, down_locs_1, up_locs_2, down_locs_2, id, sickType, chr, abpos, ref, alt, gene, extraInfo)
    end
    close(dnas)
    close(w)
end

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

main("REBASE.txt","snvs_modified.fasta","Results/step3")