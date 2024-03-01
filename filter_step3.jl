using BioSequences
using FASTX

function cal_restriction_count(dna::LongSequence{DNAAlphabet{4}}, site::LongSequence{DNAAlphabet{4}})::Int64
    query = ExactSearchQuery(site, iscompatible)
    reverse_dna = reverse_complement(dna)
    length(findall(query,dna))+length(findall(query,reverse_dna))
end

function cal_target_count(extraInfo::String, extraInfo2::String)::Int64
    ## shoule calculate how many unique sites in extraInfo2
    length(split(extraInfo,",")) + length(unique([replace(m.match,r"[+-]"=>"") for m in eachmatch(r"[+-]\d+",extraInfo2)]))
end

function filterCond(longseq::LongSequence{DNAAlphabet{4}},extraInfo::String, extraInfo2::String, site::LongSequence{DNAAlphabet{4}})::Bool
    if cal_restriction_count(longseq,site) == cal_target_count(extraInfo,extraInfo2)
        return true
    else
        return false
    end
end

function filterSiteBase(dna::String, site_name::String, site::LongSequence{DNAAlphabet{4}},outF::String)
    dnas = open(FASTA.Reader, dna, index= string(dna, ".fai"))
    targetFile = joinpath("firstStep/cond3", site_name*".tsv")
    w = open(outF, "w")
    for line in eachline(targetFile)
        id, chr, abpos, ref, alt,gene, sickType, extraInfo, extraInfo2 = split(line, '\t')
        extraInfo = convert(String,extraInfo)
        extraInfo2 = convert(String,extraInfo2)
        dna = LongSequence{DNAAlphabet{4}}(FASTA.sequence(dnas[id]))
        des = FASTA.description(dnas[id])
        pos = parse(Int64,split(des," ")[2])
        longseq = dna[(pos-120):(pos+120)]
        if filterCond(longseq,extraInfo,extraInfo2, site)
            println(w, id, '\t', chr, '\t', abpos, '\t', ref, '\t', alt, '\t', gene, '\t', sickType, '\t', extraInfo, '\t', extraInfo2)
        end
    end
    close(w)
    close(dnas)
end

function main(rebase::String,inF::String,outpath::String)
    if !ispath(outpath)
        mkpath(outpath)
    end
    for line in eachline(rebase)
        if startswith(line,"Name")
            continue
        end
        name, site = split(line, "\t")
        name = convert(String,name)
        site = LongSequence{DNAAlphabet{4}}(site)
        outname = joinpath(outpath,name * ".tsv")
        filterSiteBase(inF, name, site, outname)
    end
end

main("REs.tsv","Raw_sequences.fasta","secondStep/cond3")


