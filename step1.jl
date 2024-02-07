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


#function printResult(rawDNA::LongSequence{DNAAlphabet{4}},loc::Int64,site_length::Int64, k::String,sickType::String)
#    dna1 = rawDNA[1:(loc-1)]
#    dna2 = rawDNA[loc:(loc+site_length-1)]
#    dna3 = rawDNA[(loc+site_length):end]
#    print("$(sickType) DNA sequence in $(k) :")
#    print(dna1)
#    printstyled(dna2, bold=true, color=:orange)
#    println(dna3)
#end
#
#function printResultFull(locs::Dict{String,Vector{Int64}},locs_mutate::Dict{String,Vector{Int64}},dna::LongSequence{DNAAlphabet{4}},site_length::Int64,sickType::String,chr::String,abpos::String,ref::String,alt::String)
#    if locs != Dict("forwad"=>[0], "reverse"=>[0]) && locs_mutate == Dict("forwad"=>[0], "reverse"=>[0])
#        println("restriction site found in ", chr,":",abpos,":",ref,"=>",alt)
#        if locs["forwad"] != [0]
#            for loc in locs["forwad"]
#                printResult(dna,loc,site_length,"forwad",sickType)
#            end
#        end
#        if locs["reverse"] != [0]
#            for loc in locs["reverse"]
#                printResult(reverse_complement(dna),loc,site_length,"reverse",sickType)
#            end
#        end
#    end
#end

function printTabule(outf::IO, locs::Dict{String,Vector{Int64}},locs_mutate::Dict{String,Vector{Int64}},id::String, site_length::Int64,sickType::String,chr::String,abpos::String,ref::String,alt::String,gene::String)
    if locs != Dict("forwad"=>[0], "reverse"=>[0]) && locs_mutate == Dict("forwad"=>[0], "reverse"=>[0])
        extraInfo = ""
        if locs["forwad"] != [0]
            for loc in locs["forwad"]
                ## should calculate the absolute position
                absloc = parse(Int64,abpos)-site_length+(loc-1)
                extraInfo = extraInfo* "+" * string(absloc) * ","
            end
        end
        if locs["reverse"] != [0]
            for loc in locs["reverse"]
                absloc = parse(Int64,abpos)+site_length-(loc-1)
                extraInfo = extraInfo* "-" * string(absloc) * ","
            end
        end
        extraInfo = extraInfo[1:end-1]
        println(outf, id,"\t",chr,"\t",abpos,"\t",ref,"\t",alt,"\t",gene,"\t",sickType,"\t",extraInfo)
    end
end

function SearchRebase(dnas::String, site::LongSequence{DNAAlphabet{4}},outF::String)
    dnas = open(FASTA.Reader, dnas)
    w = open(outF, "w")
    for record in dnas
        dna = LongSequence{DNAAlphabet{4}}(FASTA.sequence(record))
        id = FASTA.identifier(record)
        id = convert(String, id)
        name = FASTA.description(record)
        chr,abpos,pos,ref,alt,_,sick,_,_,gene,_ = split(name,"|")
        chr = replace(chr,r"^\d+\s"=>"")
        sick = convert(String, sick)
        chr = convert(String, chr)
        abpos = convert(String, abpos)
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
            locs = find_restriction_site(dna, site)
            locs_mutate = find_restriction_site(mutate_dna, site)
            #printResultFull(locs,locs_mutate,dna,site_length,sick,chr,abpos,ref,alt)
            printTabule(w, locs,locs_mutate,id,site_length,sick,chr,abpos,ref,alt,gene)
        end
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
        site = LongSequence{DNAAlphabet{4}}(site)
        outname = joinpath(outpath,name * ".tsv")
        SearchRebase(inF, site, outname)
    end
end

main("REBASE.txt","snvs_modified.fasta","Results/step1")
