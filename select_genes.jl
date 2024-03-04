using DataFrames
using CSV
using DataFramesMeta
using ProgressMeter

function modifyInfo(geneF::String, ProF::String, cond::String, enzyme::AbstractString)
    # read genes
    genes=DataFrame(CSV.File(geneF,delim="\t",header=true))
    # read sample
    df=DataFrame(CSV.File(joinpath("secondStep",cond,string(enzyme,".tsv")),delim="\t",header=false))
    if ncol(df) ==8
        ## rename columns
        df= @rename df begin
            :ID = :Column1
            :chrom = :Column2
            :pos = :Column3
            :ref = :Column4
            :alt = :Column5
            :gene = :Column6
            :mutateG = :Column7
            :extrainfo = :Column8
        end
    else
        df = @rename df begin
            :ID = :Column1
            :chrom = :Column2
            :pos = :Column3
            :ref = :Column4
            :alt = :Column5
            :gene = :Column6
            :mutateG = :Column7
            :extrainfo = :Column8
            :extrainfo2 = :Column9
        end
    end
    ## read protein info
    proteins = DataFrame(CSV.File(ProF,delim="\t",header=true))
    proteins =@rename proteins begin
        :ID = :ID
        :mutateP = :mutateP
    end
    ## filter genes in sample
    df = @subset(df,in.(:gene,Ref(genes.Gene)))
    ## add protein info to df according to the ID
    df = leftjoin(df,proteins,on=:ID)
    ## add enzyme info from filename
    df = @transform(df, :enzyme = enzyme)
    ## add condition info from dirname
    df = @transform(df, :condition = cond)
    # 如果extrainfo2存在，将它放在最后一列
    if "extrainfo2" in names(df)
        df = @select(df,:gene,:mutateP,:condition,:enzyme,:ID,:chrom,:pos,:ref,:alt,:mutateG,:extrainfo,:extrainfo2)
    else
        df = @select(df,:gene,:mutateP,:condition,:enzyme,:ID,:chrom,:pos,:ref,:alt,:mutateG,:extrainfo)
    end
    df
end
function readREs(REfile::String)
    df = DataFrame(CSV.File(REfile,delim="\t",header=true))
    df.Name
end

function main(cond::String,enzyF::String,geneF::String,ProF::String)
    enzymes = readREs(enzyF)
    newL = Vector{DataFrame}(undef,length(enzymes))
    i = 1
    p = Progress(length(enzymes); dt=1, barglyphs=BarGlyphs("[=> ]"), color=:yellow)
    for enzyme in enzymes
        filepath = joinpath("secondStep",cond,string(enzyme,".tsv"))
        ## 判断文件内容是否为空
        if filesize(filepath) == 0
            newL[i] = DataFrame()
        else
            newL[i] = modifyInfo(geneF,ProF,cond,enzyme)
        end
        i += 1
        next!(p)
    end
    ## 过滤掉空的元素
    newL = filter(x->!isempty(x),newL)
    ## 按行合并所有的dataframe
    newDf = foldl(vcat, newL)
    ## 输出文件
    newDf |> CSV.write("$(cond).select.tsv",delim="\t")
end

main("cond2","REs.tsv","top_mutations.tsv","id2protein.tsv")