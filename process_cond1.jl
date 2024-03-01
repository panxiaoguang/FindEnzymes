using BioSequences
using FASTX

function find_restriction_site(dna::LongSequence{DNAAlphabet{4}}, site::LongSequence{DNAAlphabet{4}})::Dict{String, Vector{Int64}}
    query = ExactSearchQuery(site, iscompatible)
    reverse_dna = reverse_complement(dna)
    forward_check = occursin(query, dna)
    reverse_check = occursin(query, reverse_dna)
    CheckResult = Dict("forwad" => Vector{Int64}(), "reverse" => Vector{Int64}())
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



function printTabule(
    outf::IO,
    locs::Dict{String, Vector{Int64}},
    locs_mutate::Dict{String, Vector{Int64}},
    id::AbstractString,
    site_length::Int64,
    chr::AbstractString,
    abpos::Int64,
    ref::AbstractString,
    alt::AbstractString,
    gene::AbstractString,
    hgvsg::AbstractString,
)
    if locs != Dict("forwad" => [0], "reverse" => [0]) && locs_mutate == Dict("forwad" => [0], "reverse" => [0])
        extraInfo = ""
        if locs["forwad"] != [0]
            for loc in locs["forwad"]
                ## should calculate the absolute position
                absloc = abpos - site_length + 1 + (loc - 1)
                extraInfo = extraInfo * "+" * string(absloc) * ","
            end
        end
        if locs["reverse"] != [0]
            for loc in locs["reverse"]
                absloc = abpos + site_length - 1 - (loc - 1)
                extraInfo = extraInfo * "-" * string(absloc) * ","
            end
        end
        extraInfo = extraInfo[1:end-1]
        println(outf, id, "\t", chr, "\t", abpos, "\t", ref, "\t", alt, "\t", gene, "\t", hgvsg, "\t", extraInfo)
    end
end

function SearchRebase(dnas::String, mutated::String, tabfile::String, site::LongSequence{DNAAlphabet{4}}, outF::String)
    dnas = open(FASTA.Reader, dnas, index = string(dnas, ".fai"))
    mutations = open(FASTA.Reader, mutated, index = string(mutated, ".fai"))
    w = open(outF, "w")
    for line in eachline(tabfile)
        if startswith(line, "ID")
            continue
        end
        id, chr, start, _, _, wt, mut, hgvsg, sym = split(line, "\t")
        start = parse(Int64, start)
        dna = LongSequence{DNAAlphabet{4}}(FASTA.sequence(dnas[id]))
        mutate_dna = LongSequence{DNAAlphabet{4}}(FASTA.sequence(mutations[id]))
        des = FASTA.description(dnas[id])
        pos = parse(Int64, split(des, " ")[2])
        site_length = length(site)
        dna = dna[(pos-site_length+1):(pos+site_length-1)]
        mutate_dna = mutate_dna[(pos-site_length+1):(pos+site_length-1)]
        locs = find_restriction_site(dna, site)
        locs_mutate = find_restriction_site(mutate_dna, site)
        printTabule(w, locs, locs_mutate, id, site_length, chr, start, wt, mut, sym, hgvsg)
    end
    close(dnas)
    close(mutations)
    close(w)
end

function readEnzymes(rebase::String)
    enzymes = Vector{NamedTuple{(:name, :site), Tuple{String, LongSequence{DNAAlphabet{4}}}}}()
    for line in eachline(rebase)
        if startswith(line, "Name")
            continue
        end
        name, site = split(line, "\t")
        name = convert(String, name)
        site = LongSequence{DNAAlphabet{4}}(site)
        push!(enzymes, (name = name, site = site))
    end
    enzymes
end

function main(rebase::String, inF::String, mutated::String, tabfile::String, outpath::String)
    if !ispath(outpath)
        mkpath(outpath)
    end
    enzymes = readEnzymes(rebase)
    Threads.@threads for ele in enzymes
        outname = joinpath(outpath, ele.name * ".tsv")
        SearchRebase(inF, mutated, tabfile, ele.site, outname)
    end
end

main("REs.tsv", "Raw_sequences.fasta", "Mutated_sequences.fasta", "Cosmic_MutantCensus_v99_GRCh38.simple.dedup.tsv", "firstStep/cond1")


