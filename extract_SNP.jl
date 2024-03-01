function filterCond(hgvsg::AbstractString) :: Bool
    if occursin("del",hgvsg) || occursin("ins",hgvsg) || occursin("dup",hgvsg) || occursin("delins",hgvsg) || occursin(">",hgvsg)
        return true
    else
        return false
    end
end


function extract_mutations(inputF::String,outF::String)
    w = open(outF, "w")
    println(w,"ID\tCHR\tSTART\tSTOP\tSTRAND\tWT_ALLELE\tMUT_ALLELE\tHGVSG\tSYMBOL")
    ## 7/14/15/16/17/23/24/22/1
    for line in eachline(inputF)
        lineInfo = split(line, "\t")
        if filterCond(lineInfo[22])
            println(w,lineInfo[7],"\t",lineInfo[15],"\t",lineInfo[16],"\t",lineInfo[17],"\t",lineInfo[18],"\t",lineInfo[24],"\t",lineInfo[25],"\t",lineInfo[23],"\t",lineInfo[1])
        end
    end
    close(w)
end

extract_mutations("Cosmic_MutantCensus_v99_GRCh38.tsv","Cosmic_MutantCensus_v99_GRCh38.simple.tsv")
