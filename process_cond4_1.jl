using BioSequences
using FASTX

function AproximateMatch(longSeq::LongSequence{DNAAlphabet{4}}, shortSeq::LongSequence{DNAAlphabet{4}}, max_mismatches::Int64)::Dict{Int64, Vector{Int64}}
    ## 判断longSeq的长度是否为奇数，否则抛出异常
    if iseven(length(longSeq))
        throw(ArgumentError("The length of longSeq must be odd"))
    end
    ## 计算longSeq中心点的坐标位置
    center = Int64((length(longSeq) + 1) / 2)
    matches = Dict{Int64, Vector{Int64}}()
    len_long = length(longSeq)
    len_short = length(shortSeq)

    for i ∈ 1:(len_long-len_short+1)  # 滑动窗口
        mismatches = 0
        ## 记录每个mismatches距离center的距离，在前为负，在后为正
        dist = Vector{Int64}()
        for j in eachindex(shortSeq)  # 比对短序列与长序列的子序列
            if !(iscompatible(longSeq[i+j-1], shortSeq[j]))
                mismatches += 1
                push!(dist, i + j - 1 - center)
                if mismatches > max_mismatches
                    break  # 超过最大错配数，停止比对
                end
            end
        end
        if mismatches <= max_mismatches
            # 记录匹配位置
            matches[i] = dist
        end
    end
    matches
end

function FileterCond(int_list::Vector{Int64})::Bool
    # 检查是否包含0
    if 0 in int_list
        return false
    end

    # 检查所有数字的正负号是否不相同
    signs = sign.(int_list)
    if length(unique(signs)) != 1
        return false
    end

    # 绝对值列表中的最小值应该大于3
    abs_int_list = abs.(int_list) # 获取绝对值列表
    if minimum(abs_int_list) <= 3
        return false
    end

    # 如果所有检查都通过，则返回true
    return true
end

function isInLeft(int_list::Vector{Int64}) ::Bool
    signing = sign(int_list[1])
    if signing < 0
        return true
    else
        return false
    end
end

function backtrace_sequence(longseq::LongSequence{DNAAlphabet{4}}, shortseq::LongSequence{DNAAlphabet{4}}, start::Int64)::LongSequence{DNAAlphabet{4}}
    truncted_seq = longseq[start:(length(shortseq)+start-1)]
    newseq = copy(shortseq)
    for i in eachindex(truncted_seq)
        if !iscompatible(truncted_seq[i], shortseq[i])
            newseq[i] = truncted_seq[i]
        end
    end
    newseq
end

function match_sequence(seq::LongSequence{DNAAlphabet{4}}, pattern::LongSequence{DNAAlphabet{4}}, max_mismatch::Int64)
    result = Dict{String, Vector{NamedTuple{(:orient, :start, :seq), Tuple{Int64, Int64, LongSequence{DNAAlphabet{4}}}}}}()
    result["forward"] = Vector{NamedTuple{(:orient, :start, :seq), Tuple{Int64, Int64, LongSequence{DNAAlphabet{4}}}}}()
    result["reverse"] = Vector{NamedTuple{(:orient, :start, :seq), Tuple{Int64, Int64, LongSequence{DNAAlphabet{4}}}}}()
    revSeq = reverse_complement(seq)
    mt1 = AproximateMatch(seq, pattern, max_mismatch)
    for (k, v) in mt1
        if FileterCond(v)
            if isInLeft(v)
                push!(result["forward"], (orient=0, start = k, seq = backtrace_sequence(seq, pattern, k)))
            else
                push!(result["forward"], (orient=1, start = k, seq = backtrace_sequence(seq, pattern, k)))
            end
        end
    end
    if result["forward"] == []
        result["forward"] = [(orient = 0, start = 0, seq = dna"NNN")]
    end
    mt2 = AproximateMatch(revSeq, pattern, max_mismatch)
    for (k, v) in mt2
        if FileterCond(v)
            if isInLeft(v)
                push!(result["reverse"], (orient=1, start = k, seq = backtrace_sequence(revSeq, pattern, k)))
            else
                push!(result["reverse"], (orient=0, start = k, seq = backtrace_sequence(revSeq, pattern, k)))
            end
        end
    end
    if result["reverse"] == []
        result["reverse"] = [(orient = 0, start = 0, seq = dna"NNN")]
    end
    result
end

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

function get_step1_IDs(site_name::String)::Vector{String}
    first_step_file = joinpath("firstStep/cond1", site_name * ".tsv")
    ids = Vector{String}()
    for line in eachline(first_step_file)
        id, _, _, _, _, _, _, _ = split(line, '\t')
        id = convert(String, id)
        push!(ids, id)
    end
    ids
end


function printTabule(
    outf::IO,
    matched::Dict{String, Vector{NamedTuple{(:orient, :start, :seq), Tuple{Int64, Int64, LongSequence{DNAAlphabet{4}}}}}},
    mutate_dna::LongSequence{DNAAlphabet{4}},
    id::AbstractString,
    chr::AbstractString,
    abpos::Int64,
    ref::AbstractString,
    alt::AbstractString,
    gene::AbstractString,
    hgvsg::AbstractString,
    site_length::Int64,
)
    if matched != Dict("forward" => [(orient = 0, start = 0, seq = dna"NNN")], "reverse" => [(orient = 0, start = 0, seq = dna"NNN")])
        extraInfo = ""
        if matched["forward"] != [(orient = 0, start = 0, seq = dna"NNN")]
            for tp in matched["forward"]
                if find_restriction_site(mutate_dna, tp.seq) == Dict("forwad" => [0], "reverse" => [0])
                    absloc = abpos - site_length + 1 + (tp.start - 1)
                    extraInfo = extraInfo * "+" * string(absloc) * "|" * string(tp.seq) * "|" * string(tp.orient) * ","
                end
            end
        end
        if matched["reverse"] != [(orient = 0, start = 0, seq = dna"NNN")]
            for tp in matched["reverse"]
                if find_restriction_site(mutate_dna, tp.seq) == Dict("forwad" => [0], "reverse" => [0])
                    absloc = abpos + site_length - 1 - (tp.start - 1)
                    extraInfo = extraInfo * "-" * string(absloc) * "|" * string(tp.seq) * "|" * string(tp.orient) * ","
                end
            end
        end
        if extraInfo != ""
            extraInfo = extraInfo[1:end-1]
            println(outf, id, "\t", chr, "\t", abpos, "\t", ref, "\t", alt, "\t", gene, "\t", hgvsg, "\t", extraInfo)
        end
    end
end

function SearchRebase(dna::String, mutated::String, tabfile::String, site_name::String, site::LongSequence{DNAAlphabet{4}}, outF::String)
    dnas = open(FASTA.Reader, dna, index = string(dna, ".fai"))
    mutations = open(FASTA.Reader, mutated, index = string(mutated, ".fai"))
    filterID = get_step1_IDs(site_name)
    w = open(outF, "w")
    ## read the file from step1
    for line in eachline(tabfile)
        id, chr, start, _, _, wt, mut, hgvsg, sym = split(line, "\t")
        id = convert(String, id)
        if id == "ID" || (id in filterID)
            continue
        else
            start = parse(Int64, start)
            dna = LongSequence{DNAAlphabet{4}}(FASTA.sequence(dnas[id]))
            mutate_dna = LongSequence{DNAAlphabet{4}}(FASTA.sequence(mutations[id]))
            des = FASTA.description(dnas[id])
            pos = parse(Int64, split(des, " ")[2])
            site_length = length(site)
            dna = dna[(pos-site_length+1):(pos+site_length-1)]
            mutate_dna = mutate_dna[(pos-site_length+1):(pos+site_length-1)]
            matched = match_sequence(dna, site, 3)
            printTabule(w, matched, mutate_dna, id, chr, start, wt, mut, sym, hgvsg, site_length)
        end
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
        SearchRebase(inF, mutated, tabfile, ele.name, ele.site, outname)
    end
end

main("REs.tsv", "Raw_sequences.fasta", "Mutated_sequences.fasta", "Cosmic_MutantCensus_v99_GRCh38.simple.dedup.tsv", "firstStep/cond4_1")


