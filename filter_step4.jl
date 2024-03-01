function filterCond(extraInfo2::String)::Bool
    if length(unique([m.match for m in eachmatch(r"up|down",extraInfo2)])) == 1
        return true
    else
        return false
    end
end

function get_all_ids(fs::String)::Vector{String}
    ids = Vector{String}()
    for line in eachline(fs)
        id, _, _, _, _,_, _, _ = split(line, '\t')
        id = convert(String,id)
        push!(ids, id)
    end
    ids
end

function get_diff_ids(fs1::String, fs2::String)::Vector{String}
    ids1 = get_all_ids(fs1)
    ids2 = get_all_ids(fs2)
    setdiff(ids1, ids2)
end

function filterSiteBase(site_name::String,outF::String)
    targetFile1 = joinpath("firstStep/cond4_2", site_name*".tsv")
    targetFile2 = joinpath("firstStep/cond5", site_name*".tsv")
    diffIDs = get_diff_ids(targetFile1, targetFile2)
    w = open(outF, "w")
    for line in eachline(targetFile1)
        id, _, _, _, _,_, _, _ = split(line, '\t')
        id = convert(String,id)
        if id in diffIDs
            println(w, line)
        end
    end
    close(w)
end

function main(rebase::String,outpath::String)
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
        filterSiteBase(name, outname)
    end
end

main("REs.tsv","secondStep/cond4_2")

